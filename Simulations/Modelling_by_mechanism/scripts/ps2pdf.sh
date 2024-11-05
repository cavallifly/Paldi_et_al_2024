
for file in $(ls *.ps); 
do 
    echo $file

    ps2pdf $file
    
    ~/pdfcrop.pl --margin 0 ${file%.ps}.pdf _tmp_${file%.ps}.pdf #&> /dev/null
    mv _tmp_${file%.ps}.pdf ${file%.ps}.pdf
    convert -background white -alpha remove -density 600 ${file%ps}pdf ${file%ps}png
    
    rm $file #${file%ps}pdf
done
