#!/bin/bash

outdir="/Users/yaylasezginer/Documents/MATLAB/CIOOS/sensor_network/SalishSeaCast"
logfile="${outdir}/download_done.log"

mkdir -p "$outdir"
touch "$logfile"

location_names=(CentralSoG SoGeast Quadra ButeInlet Baynes)
x_array=(269 283 138 246 126)
y_array=(425 417 762 835 607)

vars=(dissolved_oxygen total_alkalinity dissolved_inorganic_carbon)

download_file () {
    location_name=$1
    X=$2
    Y=$3
    var=$4
    year=$5
    month=$6

    mm=$(printf "%02d" "$month")

    outfile="${outdir}/${location_name}_${var}_${year}_${mm}.nc"
    key="${location_name},${var},${year}-${mm}"

    if grep -q "^${key}$" "$logfile"; then
        echo "SKIP (logged): $key"
        return 0
    fi

    if [ -s "$outfile" ]; then
        echo "SKIP (exists): $outfile"
        echo "$key" >> "$logfile"
        return 0
    fi

    last_day=$(date -j -v+1m -v-1d \
        -f "%Y-%m-%d" \
        "${year}-${mm}-01" \
        +"%d")

    url="https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DChemistryFields1hV21-11.nc?${var}%5B(${year}-${mm}-01T00:30:00Z):1:(${year}-${mm}-${last_day}T23:30:00Z)%5D%5B(0.5000003):1:(441.4661)%5D%5B(${Y}):1:(${Y})%5D%5B(${X}):1:(${X})%5D"

    tmpfile="${outfile}.tmp"

    curl --fail --silent --show-error \
        --retry 6 --retry-delay 10 \
        -o "$tmpfile" "$url"

    if [ ! -s "$tmpfile" ]; then
        rm -f "$tmpfile"
        return 1
    fi

    mv "$tmpfile" "$outfile"

    echo "$key" >> "$logfile"
    echo "DONE: $outfile"
}

export -f download_file
export outdir logfile

for i in "${!location_names[@]}"; do
    location_name=${location_names[$i]}
    X=${x_array[$i]}
    Y=${y_array[$i]}

    for var in "${vars[@]}"; do
        for year in {2015..2026}; do
            for month in {1..12}; do
                echo "$location_name $X $Y $var $year $month"
            done
        done
    done
done | parallel -j 4 --colsep ' ' download_file {1} {2} {3} {4} {5} {6}
