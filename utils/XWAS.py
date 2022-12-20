# OLS or logistic regression for GWAS, TWAS or other omic-wide associations
import statsmodels.formula.api as smf
import pandas as pd

def _check_input(inputs):
    '''
    Check if input is valid:
    1. a file name (string) ending with .txt or .csv: Read in the file and return a DataFrame. Or
    2. Print error message if input is not a string
    '''
    # Read file into dataframe if not already
    if isinstance(inputs, str):
        if inputs.endswith('.txt'):
            try:
                df = pd.read_csv(inputs, sep='\t')
                return df
            except:
                print(f'#Error: File {inputs} not found\nExit')
        elif inputs.endswith('.csv'):
            try:
                df = pd.read_csv(inputs)
                return df
            except:
                print(f'#Error: File {inputs} not found\nExit')
        else:
            print('#Error: File name must ends with .txt or .csv\Exit')
    else:
        # Print error if input is not file name or dataframe
        print('#Error: Input is not a string or DataFrame\nEixt')

def _check_val_names(var_names, col_names):
    '''
    Check output name and covariate names are valid column names
    - var_name: a list of variable names
    Return False and the invalid variable name
    '''
    for val in var_names:
        if val not in col_names:
            return (False, val)
    return (True, None)

def regressor(inputs,
              outcome:str='',
              covariates:list='',
              output_fn:str='',
              reg_type:str='OLS',
              p_adj=None):
    '''
    This function performs OLS or logistic regression for GWAS, TWAS or other omic-wide associations
    Param:
    - inputs: input file name (use .txt or .csv suffix) or dataframe
    - outcome: a single column name of dependent variable (y)
    - covariates: independent variable(s) (X)
    - output_fn: output file
    - reg_type: Default is 'OLS'. OLS (continuous) or logistic regression (binary)
    - p_adj: Default is None. Multiple testing adjustment of p values. Bonferroni, FDR or None (no adjustment)
    Return:
    - a list of results: beta, p-values, SE
    '''
    # ----------------- Sanity checks -----------------
    # Make sure dependent and independent variables are present in input file/dataframe
    if not isinstance(inputs, pd.DataFrame): # Read in file if input is not a DataFrame yet
        inputs = _check_input(inputs)

    # Check if variable names exist in column names
    var_name_flag, _ = _check_val_names([outcome] + covariates, inputs.columns)
    if not var_name_flag:
        print(f'#Error: column {_} not found in input dataframe')

    # ----------------- Regression -----------------
    # Determine regression type
    # y = inputs[outcome]
    # X = inputs[covariates]
    formula = f'{outcome} ~ ' + '+'.join(covariates)

    if reg_type.upper() == 'ols':
        reg = smf.ols(formula, inputs, missing='drop') # Drop missing values when necessary
    elif reg_type == 'log':
        reg = smf.logit(formula, inputs, missing='drop')
    elif reg_type == 'mixed':
        reg = smf.mixedlm(formula, inputs, missing='drop')
    else:
        print(f'# Error: Regression type is not clear: {reg_type}')
        return
    return reg