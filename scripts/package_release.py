#!/usr/bin/env python3
"""Build/test/package a portable Linux binary; publishing is handled by CI.

Python 3.12+ standard library only. Requires Rust 1.93.1, its musl target,
musl-gcc, gcc/binutils, and the workspace's locked Cargo dependencies.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
TARGET = 'x86_64-unknown-linux-musl'
TOOLCHAIN = '1.93.1'
FLAGS = '-C target-cpu=x86-64 -C target-feature=+crt-static'


def run(args, **kwargs):
    return subprocess.run(list(map(str, args)), check=True, cwd=ROOT, **kwargs)


def output(args, **kwargs):
    return run(args, text=True, stdout=subprocess.PIPE, **kwargs).stdout.strip()


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def smoke(binary, fixtures, work):
    """Check both count policies on the shipped fixture using this exact binary."""
    index = work / 'reference.ecidx'
    run([binary, 'index-ec', '--threads', '2', '--transcripts', fixtures/'reference.fa',
         '--t2g-map', fixtures/'t2g.tsv', '--format', 'ptrhash', '--syncmer-s', '20', '--out', index])
    for mode in ('intersection', 'usa'):
        out = work / mode
        mates = [fixtures/'R1.fastq', fixtures/'R2.fastq']
        if mode == 'usa':
            mates = [work/'R1.fastq.gz', work/'R2.fastq.gz']
            for mate, source in zip(mates, (fixtures/'R1.fastq', fixtures/'R2.fastq'), strict=True):
                mate.write_bytes(gzip.compress(source.read_bytes(), mtime=0))
        run([binary, 'quant', '--threads', '2', '--index', index,
             '--r1', mates[0], '--r2', mates[1],
             '--barcode-whitelist', fixtures/'barcodes.txt', '--resolution', mode, '--out-dir', out])
        metrics = json.loads((out/'count.metrics.json').read_text())
        expected = dict(rows=7, cells=2, molecules_after_umi_correction=6,
                        spliced=4, unspliced=1, compatible_both=1)
        for key, value in expected.items():
            if metrics[key] != value:
                raise RuntimeError(f'{mode}: {key}={metrics[key]}, expected {value}')
        if metrics['correct_umis'] != (mode == 'intersection'):
            raise RuntimeError('Wrong UMI correction policy')
        if (out/'assignments.cstrad2').exists():
            raise RuntimeError('Successful quant did not clean its intermediate')


def dependency_notices(metadata, stage):
    packages = {p['id']: p for p in metadata['packages']}
    nodes = {n['id']: n for n in metadata['resolve']['nodes']}
    todo = [p['id'] for p in metadata['packages'] if p['name'] == 'constellation-cli']
    visited = set()
    while todo:
        key = todo.pop()
        if key in visited:
            continue
        visited.add(key)
        for dep in nodes[key]['deps']:
            if any(kind['kind'] != 'dev' for kind in dep['dep_kinds']):
                todo.append(dep['pkg'])
    inventory = []
    for package in sorted((packages[k] for k in visited), key=lambda p: (p['name'], p['version'])):
        if package['source'] is None:  # project license files are at archive root
            continue
        source = Path(package['manifest_path']).parent
        dest = stage/'third-party'/f"{package['name']}-{package['version']}"
        notices = []
        for path in sorted(source.rglob('*')):
            if path.is_file() and path.name.lower().startswith(('license', 'copying', 'notice', 'copyright')):
                relative = path.relative_to(source)
                target = dest/relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
                notices.append(str(target.relative_to(stage)))
        if package.get('license_file'):
            path = source/package['license_file']
            target = dest/path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            notices.append(str(target.relative_to(stage)))
        if not notices:
            curated = ROOT/'release/third-party-notices'
            entries = json.loads((curated/'provenance.json').read_text()).get(f"{package['name']}-{package['version']}", [])
            for entry in entries:
                path = curated/entry['path']
                if sha(path) != entry['sha256']:
                    raise RuntimeError(f'Upstream notice checksum mismatch: {path}')
                relative = Path(entry['path']).relative_to(f"{package['name']}-{package['version']}")
                target = dest/relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
                notices.append(str(target.relative_to(stage)))
        if not notices and package.get('license') in {'MIT', 'MIT OR Apache-2.0', 'Apache-2.0 OR MIT', 'Apache-2.0/MIT', 'MIT/Apache-2.0'}:
            # Some published crates and upstream revisions contain only the
            # SPDX declaration. Preserve their complete source/author metadata
            # and supply the declared standard text without inventing a holder.
            dest.mkdir(parents=True, exist_ok=True)
            text = (ROOT/'LICENSE-MIT').read_text().replace('Copyright (c) 2026 Constellation contributors\n\n', '')
            target = dest/'DECLARED-LICENSE-MIT.txt'
            target.write_text('This crate offers MIT in Cargo.toml. Original source and author metadata accompany this standard license text.\n\n'+text)
            notices.append(str(target.relative_to(stage)))
            shutil.copytree(source, dest/'source', ignore=shutil.ignore_patterns('.git', 'target'))
        # Include the unmodified source for MPL-covered dependencies, in addition
        # to the declared license text. Build-only dependencies are included too.
        if 'MPL' in (package.get('license') or ''):
            shutil.copytree(source, dest/'source', ignore=shutil.ignore_patterns('.git', 'target'))
        if not notices:
            raise RuntimeError(f"No license/copyright file found for {package['name']} {package['version']}")
        inventory.append({k: package.get(k) for k in ('name', 'version', 'authors', 'license', 'repository', 'source')} | {'notices': notices})
    (stage/'THIRD-PARTY.json').write_text(json.dumps(inventory, indent=2)+'\n')
    shutil.copy2(ROOT/'release/third-party-notices/provenance.json', stage/'third-party/upstream-notices.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', help='Require an exact version tag, e.g. v0.1.0 (CI release only)')
    parser.add_argument('--allow-dirty', action='store_true', help='Build a clearly labeled local development snapshot')
    parser.add_argument('--out-dir', type=Path, default=ROOT/'dist')
    args = parser.parse_args()
    version = tomllib.loads((ROOT/'crates/constellation-cli/Cargo.toml').read_text())['package']['version']
    dirty = bool(output(['git', 'status', '--porcelain', '--untracked-files=normal']))
    commit = output(['git', 'rev-parse', 'HEAD'])
    if dirty and (not args.allow_dirty or args.tag):
        parser.error('Release builds require a clean tree; --allow-dirty is only for local snapshots.')
    if args.tag:
        if args.tag != f'v{version}':
            parser.error(f'Tag must equal Cargo version v{version}')
        if output(['git', 'rev-parse', f'{args.tag}^{{commit}}']) != commit:
            parser.error('Release tag must point to the checked-out commit')
    env = dict(os.environ, RUSTFLAGS=FLAGS)
    env.pop('CARGO_ENCODED_RUSTFLAGS', None)
    # Keep the build output location unambiguous even with a user Cargo setting.
    env['CARGO_TARGET_DIR'] = str(ROOT/'target')
    run(['cargo', f'+{TOOLCHAIN}', 'build', '--locked', '--profile', 'dist', '--target', TARGET,
         '-p', 'constellation-cli'], env=env)
    binary = ROOT/'target'/TARGET/'dist/constellation'
    if 'INTERP' in output(['readelf', '-l', binary]) or '(NEEDED)' in output(['readelf', '-d', binary]):
        raise RuntimeError('Expected a static executable without a dynamic loader or shared libraries')
    if output([binary, '--version']) != f'constellation {version}':
        raise RuntimeError('Binary/Cargo version mismatch')
    metadata = json.loads(output(['cargo', f'+{TOOLCHAIN}', 'metadata', '--locked', '--format-version', '1', '--filter-platform', TARGET], env=env))
    source_files = [ROOT/'Cargo.toml', ROOT/'Cargo.lock']
    for folder in ('crates', 'examples/quickstart', 'release'):
        source_files.extend(p for p in (ROOT/folder).rglob('*') if p.is_file())
    source_files.append(Path(__file__).resolve())
    sources = {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(source_files))}
    fingerprint = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    label = f'{version}-dev-{fingerprint[:12]}' if dirty or not args.tag else version
    stem = f'constellation-{label}-{TARGET}'
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    archive = out/f'{stem}.tar.gz'
    with tempfile.TemporaryDirectory(prefix='constellation-package-') as tmp:
        work = Path(tmp)
        stage = work/stem
        stage.mkdir()
        shutil.copy2(binary, stage/'constellation')
        (stage/'constellation').chmod(0o755)
        shutil.copytree(ROOT/'examples/quickstart', stage/'examples/quickstart')
        for name in ('LICENSE-MIT', 'LICENSE-APACHE'):
            shutil.copy2(ROOT/name, stage/name)
        shutil.copy2(ROOT/'release/install.sh', stage/'install.sh')
        (stage/'install.sh').chmod(0o755)
        readme = (ROOT/'release/README.md').read_text().replace('@VERSION@', label).replace('@ARCHIVE@', archive.name).replace('@DIRECTORY@', stem)
        (stage/'README.md').write_text(readme)
        dependency_notices(metadata, stage)
        rust_docs = Path(output(['rustc', f'+{TOOLCHAIN}', '--print', 'sysroot']))/'share/doc/rust'
        shutil.copytree(rust_docs/'licenses', stage/'third-party/rust/licenses')
        for name in ('COPYRIGHT.html', 'COPYRIGHT-library.html'):
            shutil.copy2(rust_docs/name, stage/'third-party/rust'/name)
        musl_copyright = Path(os.environ.get('CONSTELLATION_MUSL_COPYRIGHT', '/usr/share/doc/musl/copyright'))
        shutil.copy2(musl_copyright, stage/'third-party/musl-copyright')
        check = work/'check'
        check.mkdir()
        smoke(stage/'constellation', stage/'examples/quickstart', check)
        info = dict(version=version, package_version=label, target=TARGET, static=True,
                    source_commit=commit, source_dirty=dirty, source_fingerprint=fingerprint,
                    source_files=sources, rustc=output(['rustc', f'+{TOOLCHAIN}', '--version']),
                    rustflags=FLAGS, binary_sha256=sha(stage/'constellation'),
                    smoke_test='passed: 7 pairs / 6 molecules for intersection and USA')
        (stage/'BUILD-INFO.json').write_text(json.dumps(info, indent=2)+'\n')
        files = sorted(p for p in stage.rglob('*') if p.is_file())
        (stage/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(stage)}\n' for p in files))
        # Normalize archive metadata; build provenance is explicit in BUILD-INFO.
        def normalize(tarinfo):
            tarinfo.uid = tarinfo.gid = 0
            tarinfo.uname = tarinfo.gname = ''
            tarinfo.mtime = 0
            return tarinfo
        with archive.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w') as tar:
                tar.add(stage, arcname=stem, filter=normalize)
        extracted = work/'extracted'
        with tarfile.open(archive) as tar:
            tar.extractall(extracted, filter='data')
        unpacked = extracted/stem
        # Test checksums and the installer from a different working directory,
        # including paths containing spaces, before considering the asset ready.
        prefix = work/'install prefix'
        run(['sh', unpacked/'install.sh', '--prefix', prefix])
        installed_check = work/'installed-check'
        installed_check.mkdir()
        smoke(prefix/'bin/constellation', unpacked/'examples/quickstart', installed_check)
    (out/f'{archive.name}.sha256').write_text(f'{sha(archive)}  {archive.name}\n')
    (out/'release-notes.md').write_text(readme)
    print(f'Ready: {archive}\nSHA-256: {sha(archive)}\nSource dirty: {dirty}')


if __name__ == '__main__':
    main()
