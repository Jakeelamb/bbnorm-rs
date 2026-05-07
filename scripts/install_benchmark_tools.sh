#!/usr/bin/env bash
set -euo pipefail

PREFIX="${PREFIX:-$HOME/.local}"
SRC_DIR="${SRC_DIR:-$PREFIX/src}"
BIN_DIR="$PREFIX/bin"
PIGZ_REPO="${PIGZ_REPO:-https://github.com/madler/pigz.git}"

mkdir -p "$SRC_DIR" "$BIN_DIR"

install_with_pacman() {
  if ! command -v pacman >/dev/null 2>&1; then
    return 1
  fi
  if command -v sudo >/dev/null 2>&1 && sudo -n true >/dev/null 2>&1; then
    sudo -n pacman -S --needed --noconfirm pigz time
    return 0
  fi
  if [[ "$(id -u)" == "0" ]]; then
    pacman -S --needed --noconfirm pigz time
    return 0
  fi
  return 1
}

install_pigz_from_source() {
  if ! command -v git >/dev/null 2>&1 || ! command -v make >/dev/null 2>&1 || ! command -v cc >/dev/null 2>&1; then
    echo "Need git, make, and a C compiler to build pigz without system package access." >&2
    return 1
  fi

  if [[ ! -d "$SRC_DIR/pigz/.git" ]]; then
    git clone --depth 1 "$PIGZ_REPO" "$SRC_DIR/pigz"
  else
    git -C "$SRC_DIR/pigz" pull --ff-only || true
  fi

  make -C "$SRC_DIR/pigz" -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
  install -m 0755 "$SRC_DIR/pigz/pigz" "$BIN_DIR/pigz"
  ln -sf "$BIN_DIR/pigz" "$BIN_DIR/unpigz"
}

if ! install_with_pacman; then
  install_pigz_from_source
fi

{
  echo "Installed benchmark tools:"
  if command -v pigz >/dev/null 2>&1; then
    printf 'pigz\t%s\t' "$(command -v pigz)"
    pigz --version
  elif [[ -x "$BIN_DIR/pigz" ]]; then
    printf 'pigz\t%s\t' "$BIN_DIR/pigz"
    "$BIN_DIR/pigz" --version
  else
    echo "pigz missing"
  fi

  if command -v unpigz >/dev/null 2>&1; then
    printf 'unpigz\t%s\t' "$(command -v unpigz)"
    unpigz --version
  elif [[ -x "$BIN_DIR/unpigz" ]]; then
    printf 'unpigz\t%s\t' "$BIN_DIR/unpigz"
    "$BIN_DIR/unpigz" --version
  else
    echo "unpigz missing"
  fi

  if [[ -x /usr/bin/time ]]; then
    printf 'gnu_time\t%s\n' /usr/bin/time
  else
    printf 'gnu_time\tmissing; scripts/measure_command.py will capture elapsed time and peak RSS\n'
  fi
}
