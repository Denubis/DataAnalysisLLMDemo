#!/usr/bin/env bash
set -euo pipefail

tar cjf - data/ | split -b 2G - original_data/data.tar.bz2.