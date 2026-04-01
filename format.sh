#!/bin/bash

for dir in include src tests;
do
    find "$dir" \
        -iname *.h -o -iname *.cpp -o -iname *.tcc \
    | xargs clang-format -style=file -i
done
