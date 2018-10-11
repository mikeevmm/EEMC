# ((?:[+-\d\.]+\s+&\s+[+-\d\.]+\s+\\\\\n?){1,20})
# \begin{supertabular}{ll}\n$0\end{supertabular}\n
pandoc -i Report.md -o Report.pdf --number-sections
if echo $* | grep -e "-v" -q
    then
        pandoc -s -i Report.md -o Report.tex --number-sections
fi