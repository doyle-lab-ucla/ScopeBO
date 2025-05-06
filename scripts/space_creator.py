import os
import pandas as pd
import numpy as np
from itertools import product as it_product
from pathlib import Path

def create_reaction_space(components, directory='./', filename='reaction_space.csv'):

    """
    Reaction scope generator
    Pass csv files with the different reaction components and their featurization. 
    Function returns a csv file with the reaction space.

    ------------------------------------------------------------------------------

    components: list
        list of csv file names (as strings). One file per starting material.
            Example: ['reactant1.csv','reactant2.csv']
        The csv files should contain the names of the compounds in the first
        column and the featurization of the compounds in the remaining columns.
        The featurization needs to be nummerical.
            Example:    name    feature1     feature2
                        A       23.1         54
                        B       5.7          80

    directory: string
        set the working directory. Default is current directory.

    filename: string
        Filename of the output csv file. Default is reaction_space.csv

    check_overwrite: Boolean
        Ask if overwritting file is ok. Default is True.
    """
    
    # Set working directory.
    wdir = Path(directory)
    # Assert that the type of components fits the requirements.
    msg="Please provide the components as a list of filenames and submit again."
    assert type(components) == list, msg

    list_of_dfs = []
    i = 0
    # Generate DataFrames for each component and add them to list_of_dfs.
    for component in components:
        csv_component = wdir.joinpath(component)

        # Assertions for existence and type of file.
        msg = "The file " + component + " was not found. Please check your input and submit again."
        assert os.path.exists(csv_component), msg
        
        # Read and clean data.
        df_component_i = pd.read_csv(csv_component, float_precision = "round_trip")
        df_component_i = df_component_i.dropna(axis='columns', how='all')

        # Set data types for all columns.
        df_component_i.iloc[:,0] = df_component_i.iloc[:,0].astype(str)
        df_component_i.iloc[:,1:] = df_component_i.iloc[:,1:].astype(float)

        # Add labels for the column names according to the reaction component.
        column_names = df_component_i.columns.tolist()
        column_names = [str(name) for name in column_names]
        column_names = ["comp" + str(i+1)+"_"+name for name in column_names] 
        df_component_i.columns = column_names
        # Add DataFrame to list.
        list_of_dfs.append(df_component_i)
        # Increase running variable for the component labelling.
        i += 1

    #Generate combinations of reactants for the reaction space.
    combinations = list(it_product(*(df.itertuples(index=False) for df in list_of_dfs)))
    # Convert to a DataFrame
    df_comb = pd.DataFrame([sum((list(row) for row in combination), []) for combination in combinations],
                      columns=sum((list(df.columns) for df in list_of_dfs), []))
    
    # Get the reactant combinations for each row.
    names_comb = []
    for _, row in df_comb.iterrows():
        # Only compound names are string values.
        names = row[row.apply(lambda x: isinstance(x, str))].tolist()
        combined_name = ""
        for i in range(len(names)):
            combined_name += names[i]
            if i+1 != len(names):
                combined_name += "."
        names_comb.append(combined_name)

    # Set the name of the reactant combination as the index. Delete columns with compound names.
    df_space = df_comb.loc[:, df_comb.apply(lambda col: not col.apply(lambda x: isinstance(x, str)).all())]
    df_space.index = names_comb

    csv_filename = wdir.joinpath(filename)  # sets name of output
    df_space.to_csv(csv_filename, index=True, mode = 'w', header=True)


    return df_space