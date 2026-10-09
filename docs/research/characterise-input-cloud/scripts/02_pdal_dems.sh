#!/usr/bin/env bash
# Full-tile ground surfaces from PDAL, at 0.25 m.
set -e
cd /tmp/opencode/cloud-research
LAS=/home/usuario/workspaces/NubePuntos/pointclouds/segmentado_RGB_color.las
mkdir -p out/dem

run () { # $1 name, $2 ... json fragment file
  echo "### $1"
  cat > out/p_$1.json <<EOF
{
  "pipeline": [
    "$LAS",
    $2,
    { "type": "writers.gdal", "filename": "/tmp/opencode/cloud-research/out/dem/$1.tif", "resolution": 0.25, "data_type": "float" }
  ]
}
EOF
  /usr/bin/time -v pdal pipeline out/p_$1.json > out/dem/$1.log 2>&1 || echo "FAILED $1"
  grep -E "Maximum resident|Elapsed" out/dem/$1.log || true
  tail -2 out/dem/$1.log
}

run smrf025 '{ "type":"filters.smrf", "cell":0.25 }' &
run csf025  '{ "type":"filters.csf", "resolution":0.25 }' &
run pmf025  '{ "type":"filters.pmf", "cell_size":0.25, "slope":1.0, "max_window_size":33 }' &
wait
ls -la out/dem/