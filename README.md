# Random geometric graphs

## hyprg -- fast embedder to hyperbolic space
### build
```
cd hyprg/
./build.sh
```

### usage
```
cd hyprg/
./hyprg --help
```

## Python utils
|script|description|
|---|---|
| analyze.py   | compute parameters of arbitrary graph |
| generate-hrg.py | generate hyperbolic random graph |
| eval-embedding.py | estimate quality of graph embedding |
| fit.py | compute graph embedding to hyperbolic space (deprecated, use ./hyprg/hyprg instead)|

### Required packages
python-networkx
python-matplotlib
python-tk
python-scipy
python-sklearn

## Python libraries (deprecated)
|library|description|
|---|---|
|graph.py|TODO refactor|
|embedding_models.py|TODO|
|loss_functions.py|TODO|
|optimization.py|TODO|
|pair_generators.py|generate labelled pairs of vertices|

