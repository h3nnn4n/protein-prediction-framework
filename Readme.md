# PPF

PPF (Protein Predition Framework) is a framework based on Rosetta for researching
ab inition protein tertiary structure predction.

# How to use:

- `git clone https://github.com/h3nnn4n/protein-prediction-framework`
- `cd protein-prediction-framework`
- `git submodule update --init`
- `cd src/external`
- `git submodule update --init`
- `cd tmscore/src`
- `gfortran TMscore.f -o TMscore`
- `cd ../../../`
- `cd src/de`
- `gfortran spicker.f -o spicker`
- `python -m pip install mock`
- `python -m pip install multiprocessing`
- `python -m pip install numpy`
- `python -m pip install pytest`
- `python -m pip install types`
- `python -m pip install uuid`
- `python -m pip install yaml`
- `python -m pip install scipy`
- `python3 main.py`

# LICENSE

All files in this repository are released under GPLv3 license, unless noted otherwise.
See [LICENSE](LICENSE) for more details.

