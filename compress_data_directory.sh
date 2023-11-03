#!/usr/bin/env bash
set -euo pipefail
rm original_data/data.tar.bz2.??
tar cjf - data/ | split -b 2G - original_data/data.tar.bz2.