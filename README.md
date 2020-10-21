# AC-CFM
AC cascading failure model based on MATPOWER for resilience analysis of power networks.

<!-- GETTING STARTED -->
## Getting Started

The following steps will guide you through getting AC-CFM run on your computer.

### Prerequisites

* Matlab 2019a or later is recommended to run AC-CFM.

* Matpower 7.1 or later is required to run AC-CFM. Please follow the [instructions on the Matpower Website for installation](https://matpower.org/about/get-started/).

* IPOPT is recommended to install, as it tends to have higher convergence for AC optimal power flows. Instructions on how to install IPOPT can be found [on the Matpower website](https://matpower.org/download/optional-solvers/). The [OPTI Toolbox](https://www.inverseproblem.co.nz/OPTI/) is a simple way of installing IPOPT.

* Test the installation of Matpower and IPOPT by running the following command in Matlab:
```
runopf(case9, mpoption('opf.ac.solver', 'IPOPT'))
```

### Installation

1. Clone the repository
```sh
git clone https://github.com/mnoebels/AC-CFM.git
```
2. Add the AC-CFM folder to your Matlab path setting (menu Home -> Set Path)



<!-- USAGE EXAMPLES -->
## Usage

Use the following example to model the cascade following a single line outage in the IEEE 39-bus test network:
```
% load default settings
settings = get_default_settings();
 
% enable verbose output – just for testing
settings.verbose = 1;
 
% model outage of line 9
initial_contingency = 9;
 
% apply the model
result = accfm(case39, struct('branches', initial_contingency), settings);
```
This should produce the following output:
```
  Demand increased by 0.1% (limit is 15.0%) and generation capacity is met. Distribute slack generation.
  Exceeded line ratings: 6-11

  Demand increased by 0.5% (limit is 15.0%) and generation capacity is met. Distribute slack generation.
  Q outside limits at generators at buses 34
  Exceeded line ratings: 3-4 3-18 10-13 13-14 14-15 16-17 17-18

   3 islands and 2 isolated nodes detected
   Island: [ 1 2 3 4 5 6 7 8 9 17 25 26 27 28 29 30 31 37 38 39 ]
    Demand increased by 20.9% (limit is 15.0%) or generation capacity is not met. Perform underfrequency load shedding of 12.2%.

    Q outside limits at generators at buses 31 39
    Voltage outside limits at buses 4 7 8
    Undervoltage load shedding applied at buses 4 7 8
    Exceeded line ratings: 1-2 8-9 9-39

     3 islands and 1 isolated nodes detected
     Island: [ 2 3 17 25 26 27 28 29 30 37 38 ]
      Demand decreased by 42.6% (limit is 15.0%). Tripping 2 smallest generators.

      Loads shed (22.54%) due to voltage collapse at buses 26 29

      Demand increased by 0.0% (limit is 15.0%) and generation capacity is met. Distribute slack generation.
      Exceeded line ratings: 2-25 2-30

        2 islands and 1 isolated nodes detected
        Island: [ 17 25 26 27 28 29 37 38 ]
         No generation available.

        Island: [ 2 3 ]
         No generation available.

        Island: [ 30 ]
         Demand decreased by 100.0% (limit is 15.0%). Tripping 1 smallest generators.

     Island: [ 4 5 6 7 8 31 ]
      Demand increased by 400.4% (limit is 15.0%) or generation capacity is not met. Perform underfrequency load shedding of 79.0%.


     Island: [ 1 39 ]
      Demand decreased by 4.0% (limit is 15.0%). Distribute slack generation.

     Island: [ 9 ]
      No generation available.

   Island: [ 15 16 19 20 21 22 23 24 33 34 35 36 ]
    Demand decreased by 8.4% (limit is 15.0%). Distribute slack generation.

   Island: [ 10 11 12 13 32 ]
    Demand decreased by 98.7% (limit is 15.0%). Tripping 1 smallest generators.

    No generation available.

   Island: [ 14 ]
    No generation available.

   Island: [ 18 ]
    No generation available.

Cascade halted. Elapsed time: 6.51s
Total load shedding: 45.05%
Load shedding UFLS: 20.97% 
Load shedding UVLS: 0.88% 
Load shedding VCLS: 4.61% 
Load shedding non-converging OPF: 0.00% 
Load shedding tripped: 18.59%
```


<!-- Troubleshooting -->
## Troubleshooting

If you run into problems, for instance exceptionally large load shedding or large amount of OPF load shedding, try the following:

1. Set PMIN (column 10) of the gen matrix to zero.
2. Increase the reactive power limits (columns 4 and 5) of the gen matrix (e.g. to +80% and -40% of PMAX).
3. Loads should be fixed (not dispatchable) before passed to AC-CFM. You can use the disp2load function provided to convert dispatchable loads to fixed loads.
