# fbWecCntrl

**NOTE: this code is fully-functional, but not supported and users should not expect responses to issues**

fbWecCntrl is a set of MATLAB functions and scripts demonstrating a causal impedance matching approach to wave energy converter (WEC) control design.
The methods applied in this code are detailed in the following paper, and is a fork of code originally published on [MHK-DR](https://mhkdr.openei.org/submissions/315).

```bibtex
@Article{Coe2020a,
  author       = {Ryan G. Coe and Giorgio Bacelli and Dominic Forbush},
  title        = {A practical approach to wave energy modeling and control},
  date         = {2020, submitted},
  journaltitle = {Renewable and Sustainable Energy Reviews},
}
```

## Getting started
This code has been tested on MATLAB 2020a (9.8.0.1323502) and has the following dependencies.

Dependency                          | Website                                                         	 | Required?
----------------------------------- | ------------------------------------------------------------------ | ---------
MATLAB                              | https://www.mathworks.com/products/matlab.html                  	 | yes
MATLAB Control System Toolbox 		  | https://www.mathworks.com/products/control.html 					         | yes
MATLAB Optimization Toolbox         | https://www.mathworks.com/products/optimization.html            	 | yes
WAFO<sup>1</sup>                    | https://github.com/wafo-project/wafo                            	 | no
export_fig<sup>2</sup>   			      | https://github.com/altmany/export_fig 							               | no

<sup>1</sup>_[WAFO](https://github.com/wafo-project/wafo) is used to produce wave spectra._ 
<sup>2</sup>_[export_fig](https://github.com/altmany/export_fig) is used to produce PDFs figures._ 

1. **Download the fbWecCntrl software**: Clone the repository or download an archive. If required, unzip the archive to a path of your choosing.

2. **Add fbWecCntrl to your MATLAB path**: Add the necessary directories to your MATLAB path using the MATLAB command prompt.

    ```matlab
    >> addpath(genpath('/path/to/fbWecCntrl'));
    >> savepath;
   ```

3. **Run the demo cases**: Run the demo cases for the WaveBot and FOSWEC.

    ```matlab
    >> demo_waveBot
    >> demo_foswec
   ```
## Software License

Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
 
fbWecCntrl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

fbWecCntrl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with fbWecCntrl.  If not, see <https://www.gnu.org/licenses/>.
