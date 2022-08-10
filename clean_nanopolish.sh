echo "## Nanopolish cleaning\n"
zcat $1 | grep -f "${2}_readnames.temp" > "${2}_greped.txt"
no_recog = $(grep -v -f "${2}_greped.txt" "${2}_readnames.temp" | uniq | wc -l)
tot = $(wc -l "${2}_readnames.temp")
echo "Total read_name number: ${tot}"
echo "Not recognized read_name number: ${no_recog}"
cut -f12,13,14,15,16,17,18,19,20 "${2}_greped.txt" > "${2}_extracols.txt"
sed -i "/^[[:space:]]*\$/d" "${2}_extracols.txt"
grep -v -f "${2}_extracols.txt" "${2}_greped.txt" > "${2}_new.txt"
mv "${2}_new.txt" "${2}_greped.txt"
n_extracol = $(grep . "${2}_extracols.txt" | wc -l )
echo "number of row with extra col: ${n_extracol}"
