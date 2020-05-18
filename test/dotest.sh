#!/bin/bash

echo "Launching test, this may take a while, processing 25k molecules..." 

time ../Lilly_Medchem_Rules.rb example_molecules.smi > okmedchem.smi

if [ ! -s okmedchem.smi ]
then
  echo "Invocation failed, no output, check build" >&2
  exit 1
fi

failures=0

for stem in bad0 bad1 bad2 bad3 okmedchem
do
  computed="${stem}.smi"
  correct="${stem}.correct.smi"

  if [ ! -s "$correct" ]
  then
    echo "Correct file '${correct} missing or empty, incomplete package" >&2
    exit 2
  fi

  if [ ! -s "$computed" ]
  then
    echo "Computation failed, did not produce '${computed}'" >&2
    exit 3
  fi

  diff -w $correct $computed

  if [ $? -ne 0 ]
  then
    echo "Failure on '${stem}'" >&2
    let failures++
  fi
done

if [ $failures -gt 0 ]
then
  echo "${failures} failed tests" >&2
else
  echo "All tests successful" >&2
  rm okmedchem.smi
  rm ok?.log
  rm bad?.smi
fi
