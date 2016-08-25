#!/bin/bash

if [ $# -eq 1 ]; then

    filebase="$1"

    echo "### pdfseparate -f 1 -l 1 ${filebase}.pdf tmp.pdf"
    pdfseparate -f 1 -l 1 ${filebase}.pdf tmp.pdf
    echo "### pdfcrop tmp.pdf ${filebase}_crop.pdf"
    pdfcrop tmp.pdf ${filebase}_crop.pdf
    echo "### rm -f tmp.pdf"
    rm -f tmp.pdf

else

    echo "### Please input a filebase"

fi
