# **Weighted Birkhoff Averages Accelerate Data-Driven Methods**

This repository contains MATLAB scripts to reproduce the data and figures from [Weighted Birkhoff Averages Accelerate Data-Driven Methods](https://arxiv.org/abs/2504.00729) by Maria Bou-Sakr-El-Tayar, [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/), and [Matthew Colbrook](https://www.damtp.cam.ac.uk/user/mjc249/home.html).

## **Paper Abstract**
Many data-driven algorithms in dynamical systems are built on ergodic averages that converge slowly. A simple idea changes this: taper the ends of the sum. Weighted Birkhoff averages can converge far faster -sometimes superpolynomially or even exponentially- and can be built directly into existing methods. We illustrate this by introducing five weighted algorithms: weighted Dynamic Mode Decomposition (wtDMD), weighted Extended DMD (wtEDMD), weighted Sparse Identification of Nonlinear Dynamics (wtSINDy), weighted spectral measure estimation, and weighted diffusion forecasting. Across examples ranging from fluid flows to El Nino data, the message is clear: weighting costs nothing, is easy to implement, and often gives better results from the same data.

## **Required Data**
Many of the data files are too large to put on GitHub. Therefore, you can download all necessary data from google drive: [https://drive.google.com/drive/folders/1mz0lm0a49nZQ8MsghGh7-qHAc45uLu6V?usp=sharing](https://drive.google.com/drive/folders/1mz0lm0a49nZQ8MsghGh7-qHAc45uLu6V?usp=sharing)

## **Repository Contents**
This repository currently contains five folders, each associated to a subsection of the demonstrations section in the paper. They are organized as follows:

- **Section 3.1**: [**Dynamic Mode Decomposition**](https://github.com/jbramburger/weighted_methods/tree/main/Dynamic%20Mode%20Decomposition). This folder contains the code to implement weighted dynamic mode decomposition (wtDMD) and compare its performance to standard dynamic mode decomposition. The method is applied to fluid flow around a cylinder data and the script can be used to load and recreate all data from the paper.
- **Section 3.2**: [**Extended Dynamic Mode Decomposition**](https://github.com/jbramburger/weighted_methods/tree/main/Extended%20Dynamic%20Mode%20Decomposition). This folder contains the code to implement weighted extended dynamic mode decomposition (wtEDMD) and compare its performance to standard extended dynamic mode decomposition. The method is applied to synthetic data gathered from the standard map on the torus and the script can be used to load and recreate all data from the paper.
- **Section 3.3**: [**Model Identification**](https://github.com/jbramburger/weighted_methods/tree/main/Model%20Identification). This folder contains code to simulate solitons in the nonlinear Schrodinger equation with a parabolic trap and extract their centre of mass. Then, we use this centre of mass data to identify a reduced-order model of the soliton dynamics using model identification. The model identification script compares EDMD, wtEDMD, the sparse identification of nonlinear dynamics (SINDy) method, and weighted SINDy (wtSINDy).
- **Section 3.4**: [**Spectral Measures**](https://github.com/jbramburger/weighted_methods/tree/main/Spectral%20Measures). This folder contains code to approximate spectral measures associated to the Koopman operator for measure-preserving systems from data. Density reconstruction uses the Chebfun package which can be download from [here](https://www.chebfun.org/download/).
- **Section 3.5**: [**Diffusion Forecasts**](https://github.com/jbramburger/weighted_methods/tree/main/Model%20Identification).

 
