#!/bin/sh


ORDER_MAX=7

if [[ $# -eq 1 ]]; then
  ORDER_MAX=$1
fi  

for i in $(seq 2 ${ORDER_MAX}); do
  echo $i
  ./hermite $i 0 0 1 "ho${i}" > /dev/null
  nl -v $i -i 0 "ho${i}_x.txt" | paste - "ho${i}_w.txt" > "ho${i}_xw.txt"
done

# row-bind
echo "order\tabscissa\tweight" | cat - ho[0-9]*_xw.txt >  hoAll_xw.txt

rm ho[0-9]*_*txt

echo "~Fine~"

