# **Weighted Birkhoff Averages Accelerate Data-Driven Methods**

This repository contains MATLAB scripts to reproduce the data and figures from [Weighted Birkhoff Averages Accelerate Data-Driven Methods](https://arxiv.org/abs/2504.00729) by Maria Bou-Sakr-El-Tayar, [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/), and [Matthew Colbrook](https://www.damtp.cam.ac.uk/user/mjc249/home.html).

## **Paper Abstract**
Many data-driven algorithms in dynamical systems are built on ergodic averages that converge slowly. A simple idea changes this: taper the ends of the sum. Weighted Birkhoff averages can converge far faster -sometimes superpolynomially or even exponentially- and can be built directly into existing methods. We illustrate this by introducing five weighted algorithms: weighted Dynamic Mode Decomposition (wtDMD), weighted Extended DMD (wtEDMD), weighted Sparse Identification of Nonlinear Dynamics (wtSINDy), weighted spectral measure estimation, and weighted diffusion forecasting. Across examples ranging from fluid flows to El Nino data, the message is clear: weighting costs nothing, is easy to implement, and often gives better results from the same data.

## **Required Data**
Many of the data files are too large to put on GitHub. Therefore, you can download all necessary data from google drive: [https://drive.google.com/drive/folders/1mz0lm0a49nZQ8MsghGh7-qHAc45uLu6V?usp=sharing](https://drive.google.com/drive/folders/1mz0lm0a49nZQ8MsghGh7-qHAc45uLu6V?usp=sharing)

## **Repository Contents**
This repository currently contains five folders, each associated to a subsection of the demonstrations section in the paper. They are organized as follows:

- **Section 3.1**: [**Dynamic Mode Decomposition**](https://github.com/jbramburger/weighted_methods/tree/main/Dynamic%20Mode%20Decomposition). This folder contains the code to implement weighted dynamic mode decomposition (DMD) and compare its performance to standard DMD. The method is applied to fluid flow around a cylinder data and the script can be used to load and recreate all data from the paper.
- **Section 3.2**: [**Extended Dynamic Mode Decomposition**](https://github.com/jbramburger/weighted_methods/tree/main/Dynamic%20Mode%20Decomposition).

 
