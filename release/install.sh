#!/bin/sh
set -eu

prefix="${HOME:?HOME is required}/.local"
case "${1:-}" in
  --prefix)
    [ "$#" -eq 2 ] && [ -n "$2" ] || { echo 'Usage: sh install.sh [--prefix DIRECTORY]' >&2; exit 2; }
    prefix=$2 ;;
  --help|-h)
    echo 'Usage: sh install.sh [--prefix DIRECTORY]'
    echo 'Installs bin/constellation under DIRECTORY (default: $HOME/.local).'
    exit 0 ;;
  '') [ "$#" -eq 0 ] || exit 2 ;;
  *) echo 'Usage: sh install.sh [--prefix DIRECTORY]' >&2; exit 2 ;;
esac

if [ "$(uname -s)" != Linux ] || [ "$(uname -m)" != x86_64 ]; then
  echo 'This archive is for Linux x86_64. Use a matching release or build from source.' >&2
  exit 1
fi
package_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
(
  cd "$package_dir"
  sha256sum --check --quiet MANIFEST.sha256
)
mkdir -p -- "$prefix/bin"
install -m 755 -- "$package_dir/constellation" "$prefix/bin/constellation"
"$prefix/bin/constellation" --version
printf 'Installed to %s/bin/constellation\nEnsure %s/bin is on PATH.\n' "$prefix" "$prefix"
