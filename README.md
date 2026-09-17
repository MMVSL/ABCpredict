# ABCpredict

ABCpredict is an AI-based tool for predicting three pharmacokinetic endpoints: human intestinal absorption (expressed as Caco-2 permeability), blood-brain barrier permeability and plasma protein binding (expressed as fraction unbound).
This repository provides the necessary Python scripts to perform predictions using pretrained machine learning models and a predefined computational environment.

## Citation
If you use ABCpredict, please cite:

*Citation information will be available soon.*

## Installation

To set up the required environment, use the provided environment file:

```sh
conda env create -f env_ABCpredict.yml
```

Then, activate the environment:

```sh
conda activate ABCpredict
```

Finally, download the ZIP archive containing the trained models:
```sh
wget http://www.mmvsl.it/wp/wp-content/uploads/PUT_LINK_HERE
```
Alternatively, you can download the zipped folder by simply opening the link http://www.mmvsl.it/wp/wp-content/uploads/PUT_LINK_HERE in your browser and saving it manually.

After downloading, extract the contents of the archive:
```sh
unzip models.zip
```

Make sure that the folder `models` is located in the repository root directory (i.e., the same folder where `ABCpredict.py` is located).

- If you used `wget` and `unzip` from within the repository root directory, the folder should already be in the correct location.
- Otherwise, move it manually into the repository root.

  For example:
  ```sh
  mv ~/Downloads/models path/to/this/repo/
  ```

After setup, your folder should look like this:
```sh
ABCpredict/
├── ABCpredict.py
├── env_ABCpredict.yml
├── example_input.csv
├── models      <--- extracted folder containing trained models
│   ├── caco2_GP_hybrid.dump
│   ├── caco2_scaler.dump
│   ├── classBBB_RF_PubChem-FP.dump
│   ├── fu_GP_RDKit-desc.dump
│   └── fu_scaler.dump
├── models.zip  <--- downloaded zipped folder
├── PubChemFingerprints.py
├── rdkit_208-desc_list.dump
└── README.md

```

## Usage

Run the predictions with Python, specifying the required inputs, for example:

```sh
python ABCpredict.py -in example_input.csv -e caco2 -out output.csv
```
Available endpoints are: 
- `caco2` for Caco-2 permeability
- `bbb` for BBB permeability
- `fu` for fraction unbound

### Input
- The input file must be in CSV format and contain a column named `SMILES` containing the molecular structures. 
- An example input CSV file (`example_input.csv`) is included in the repository for testing purposes.

### Output
- The ABCpredict.py script will generate a CSV file containing predicted value for the selected endpoint.
