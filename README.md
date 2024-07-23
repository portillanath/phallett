# phallett

🎨 Phallett is a Viral Genomic Distance Profiler. Unified genomic distances metrics calculation and create graphs to visualize patterns and correlation following the guidelines of the ICTV for viral taxonomy throught the access to de Virus Metada Resouce VMR 👾. So, you can used to send your next Taxa proposoal (TaxoProps)!

## Installation. 🌈

* Linux/Ubuntu
1. You must have anaconda or miniconda, if that is not the case follow the next link: 
https://shorturl.at/ZDEL0

2. Clone this repo: `git clone https://github.com/portillanath/phallett.git`

3. Move inside the repo folder: `cd phallett`

4. Set up depedencies throught conda enviroment `bash enviroments.sh`

## It is a executive wrap 🌯
The structure of phallett are based on modules, each module is a exacutable work of tool or analysis. By default all the workflow is execute if there is not sumministed any module name `bash phallet.sh`

Modules here:

| Flag          | Action                          |
| ------------- | ------------------------------- |
|    `-d`       |  Data used with `taxa` module, default test_genus.txt  |
|    `-m`       | Module selection                |
|   `-m ictv`   | ICTV module class display the VRM update|
|   `-m taxa`   |  Selection of taxa file to analize |
|   `-m mash`   | Compute mash distances with mash and sourmash algorithms |
| `-g`          | Genus to calculate metric |
| `-km`         | Kmers for mash metrics, separate with space|
| `-m ani`      | Compute ani metrics with fastani and skani algorithms|
| `-ka`         | Kmers for ani metrics separate with space  |
| `-f`          | Fragment length for ANI default is 500     |
| `-my`         | Metric in y scatter axis, default is mash  |
| `-mx`         | Metric in x scatter axis, deafult is ani   |
| `-kx`         | Kmers used on x metric, separate with comma |
| `-ky`         | Kmers used on y metric, separate with comma|
| `-m wraggling`| Arrange summary metrics to csv             | 
| `-m graphs`    | Sccatter graphs are exported as PDF       |
