echo "# Workflow Run info:"
echo $(nextflow log | grep "\s${1}\s")
nextflow log $1 -t template_nextflowlog.md > log_temp.txt
echo "# Processes info:"
grep -o 'Arguments.*' log_temp.txt | head -n1
echo "# Processes outputs:"
grep -v 'Arguments.*' log_temp.txt
rm -f log_temp.txt
