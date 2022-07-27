#!/bin/bash

var="$(grep -n size_SC= RunNextGen.jl | cut -d : -f 1)"
echo $var

new_number=20
sed -i '$var s/.*/size_SC=10/' RunNextGen.jl
