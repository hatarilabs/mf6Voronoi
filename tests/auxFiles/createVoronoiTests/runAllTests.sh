#!/bin/bash
echo "Removing and creating output folder"
rm -rf "output"
mkdir "output"

# Run Python script 1
echo "Running test 1"
python3 testEx1RegionalModel.py

# Run Python script 2
echo "Running test 2"
python3 testEx2RiverAquifer.py

# Run Python script 3
echo "Running test 3"
python3 testEx3SiteDewatering.py

# Run Python script 3
echo "Running test 4"
python3 testEx4TrenchExcavation.py

# Run Python script 3
echo "Running test 5"
python3 testEx5stibniteMine.py
# You can add more scripts as needed
