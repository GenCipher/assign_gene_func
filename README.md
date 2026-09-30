# Part 1: Decoding gene function
When annotating candidate genes or ORFs on a new genome, some may correspond to actual genes and some detected by chance. How can we tell which ORFs correspond to the actual genes? And how can we tell what the proteins are and do? 

**The tasks of this assignment are:**
- Implement Global alignment (Needleman-Wunsch)
- Choosing our reference genomes 
- Implement Local alignment (Smith-Waterman) 
- Finding homologous genes

Open `gene_func_workbook.ipynb` to get started.

## Submission and Grading
1. *Coding exercises* require you to implement an algorithm. We will always include Python stubs for the functions you need to implement in a separate file (most often `helper_functions.py`). Please read the function docstrings for expected parameter and return types. Coding problems will be automatically graded with unit tests.

2. *Image answers* require you to generate an image and save the plot to a corresponding file in `<img>.png`. Please ensure that your images have the exact same name as specified in the instructions and are stored in the root folder of the repository (the same folder as the README.md file).

3. *Variable answers* require you to write down your answers into variables. Please be careful that the variables have the exact same name as in the instructions. We will always provide a stub in the notebook. We don't actually run your entire notebooks, we evaluate only the variable we are grading. This means that your variables should be set explicitly. For instance `x = 5 / 3` or `answer_var = x` **will not work** and will receive zero points. All answer variables should be set explicitly e.g. `answer_var = 1.6666` in for this simple example. Please ensure there are no syntax errors in your notebook. We've added GitHub Actions that test syntax errors when you push your code to GitHub. If there is a green checkmark at the top of your repository, your code has no syntax errors. If there is a red X at the top, your notebook did not compile correctly, so make sure you eliminate any errors.

## Environment instructions

You will need Python 3.10 or higher. You will need to install `biopython` for accessing NCBI and `matplotlib` for plotting. You will also need `jupyterlab` to open and run the notebook. You can also use `numpy`, `pandas`, `seaborn`, and `tqdm`. You can install everything necessary by running
```
pip install biopython matplotlib jupyterlab
```
Do not use any other libraries at this point.

You can start the notebook by running
```
jupyter lab
```
and a browser window should pop open, or you can manually navigate to http://localhost:8888/ in your browser.

**Note, this has been adapted from homework-2 from the IB-ULFRI course.** 
[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
