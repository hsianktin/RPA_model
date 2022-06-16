# Dependencies
This repository is dependent fully on julialang and its ecosystem.
Modules used
- [Statistics](https://github.com/JuliaLang/Statistics.jl)
- Random
- [CSV](https://github.com/JuliaData/CSV.jl)
- [Plots](https://doi.org/10.5281/zenodo.4725317)
- [DataFrames](https://github.com/JuliaData/DataFrames.jl)
- Printf
- [ProgressMeter](https://github.com/timholy/ProgressMeter.jl)
- Distributed
- DelimitedFiles

# Initialization
Make sure you have julialang and relevant packages installed.
Clone this package, and run `julia initialize.jl` in the command line. This is to precompile the package and generate the folders.

# Experiment Data Preparation
The functions for reading experiment data are written in [[evaluate_base]] and [[evaluate_base_diffusion]]. 

As is mentioned in the paper, we prepare the data by first incubating ssDNA with the RPA for 30 minutes and changed to a higher concentration of RPA later. This is the situation our code is primarily intended to analyze.

Data are stored in stored in the form of csv files. Each csv file consists of a single column of DNA extensions recorded by the ssDNA curtain assay.  The first row is ssDNA extension measured at time=5 seconds after flow. The n-th row is ssDNA extension measured at time $5n$ seconds after flow. 

```csv
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
1.000000000000000000e+00
...
```

The data are first grouped by salt, protein concentrations and type of proteins. Each elementary group of data are placed under the same folder. Elementary groups are further grouped together by salt concentration and type of proteins. The secondary groups are a basic set of experiments

The directory tree for the experiment data folder is given as follows.
```
.
├── data_exp
│   ├── exp_data_lwt0
│   ├── exp_data_lwt1
│   ├── exp_data_lwt10
│   ├── exp_data_lwt25
│   ├── exp_data_lwt4
│   ├── exp_data_lwt50
│   ├── exp_data_wt0
│   ├── exp_data_wt1
│   ├── exp_data_wt10
│   ├── exp_data_wt25
│   ├── exp_data_wt4
│   └── exp_data_wt50

```

The mapping from folder name to experiment conditions (experiment labels) is stored as a dictionary in the [[evaluate_base]].
```julia
    if exp_label == "wt_15mM_salt"
        global exp_data_base=[("wt",0),("wt",1),("wt",4),("wt",10),("wt",25),("wt",50)]
        global folds = [0,1,4,10,25]
    ...
    end
    function salt_concentration(type) 
        if type == "wt"
            return "15mM salt"
        elseif type == "lwt"
            return "150mM salt"
        end
    end
    function exp_dict_inject(type::String,concentration::Float64)
        datapath = "$(mainpath)/data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$concentration"]=[datapath,fig_label]
    end
    function exp_dict_inject(type::String,concentration::Int)
        datapath = "$(mainpath)/data_exp/exp_data_$(type)$(concentration)/"
        fig_label = "$(salt_concentration(type)) yRPA 0.1nM to $(concentration)nM"
        exp_dict["$(convert(Float64,concentration))"]=[datapath,fig_label]
    end
    function exp_dict_inject(tuple::Tuple)
        if length(tuple) == 2
            type,concentration = tuple
            exp_dict_inject(type,concentration)
        else
            println("input not match")
        end
    end
    exp_dict_inject.(exp_data_base)
```

# Topology of Codes and Datafiles
## dependency
- [[TonksGaswithReactionGaps]] -> [[simu_base]] -> [[auto_update]]/[[sensitivity_length]]/[[sensitivity_diffusion]]/[[sensitivity_parameter]]
- [[evaluate_base]] -> [[evaluate]]/[[evaluate_sensitivity]]

# Basic Usage
At the root directory, simulation data can be obtained by executing the following command in powershell or bash:
```Shell
julia simu_base.jl $k_on $k_off $v_open $v_close $fold $L $T1 $T2 $N $exp_label $simu_label $gaps_type
```
where `$k_on`, `$k_off`, ..., should be replaced by numerical values. `exp_label` and `simu_label` are strings and concatenated to create a unique identifier for the output. `$gaps_type` can be `exact`, `cumulative` or `none`. It corresponds to different ways of couting the gaps, or not counting them at all.

All the outputs are stored in the folder `data_simu/`. To generate trace figures, use the following command:
```shell
julia evaluate.jl $exp_label $simu_label
```
This command will combine the simulation and experimentally obtained traces. The code in [[evaluate]] could be tailored to custom needs. Accompanying the figures, the code also generates the data source for the figures.

To analyze the gap distribution, use [[exact_gap_analysis]] as follows:
```Shell
julia exact_gap_analysis.jl $exp_label $simu_label
```

# Batch Analysis
[[auto_update]] is built on top of the basic analysis method, to automatically find an optimal value of parameters. It relies on the file `figs/sources/para.csv` to get the inital parameters.

[[auto_summarize]] summarizes the data generated by [[auto_update]] and the output parameters will be saved in `figs/para_$(simu_label).csv`
The default `simu_label` for [[auto_summarize]] is `fitted`.

[[sensitivity_parameter]], [[sensitivity_length]], [[sensitivity_diffusion]] detect the local sensitivity of the model at the optimal data points. In those cases, only one parameter is varied and we inspect the losses as a function of this parameter.