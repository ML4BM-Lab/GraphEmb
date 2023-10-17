import helper_functions as hf
import pandas as pd
from upsetplot import generate_counts
from upsetplot import plot, UpSet
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': hf.get_dtis_drugbank(),
             'Davis_et_al': hf.get_dtis_davis(),
             'BindingDB': hf.get_bindingdb_dtis(),
             'BIOSNAP': hf.get_biosnap_dtis()
             }

dict_dfs.update(hf.get_dict_dtis_yamanishi())

def plot_upset(mode):

    def clean_plots():
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()

    if mode == 'drugs':
        column = 0
    else:
        column = 1

    ###
    datasets = ['DrugBank','BindingDB','BIOSNAP','Davis_et_al','E','IC','GPCR','NR']

    ## generate upset
    element_dict = {}
    for data in datasets:
        element_dict[data] = set(dict_dfs[data].iloc[:,column])

    # Compute pairwise intersections
    intersections = {}
    keys = list(element_dict.keys())
    for i in range(len(keys)):
        intersections[keys[i]] = len(element_dict[keys[i]].difference(set.union(*[element_dict[k] for k in keys if k != keys[i]])))
        for j in range(i + 1, len(keys)):
            key1, key2 = keys[i], keys[j]
            intersection = (element_dict[key1].intersection(element_dict[key2])).difference(set.union(*[element_dict[k] for k in keys if k not in (key1,key2)]))
            intersections[f'{key1} ∩ {key2}'] = len(intersection)

    # Create a DataFrame from the intersection data
    df = pd.DataFrame(intersections, columns = intersections.keys(), index = [0]).T

    def return_boolean_arrays(key, datasets):
        return tuple([True if dt in key else False for dt in datasets])

    # Create a multi-index DataFrame
    df.index = pd.MultiIndex.from_tuples([return_boolean_arrays(key, datasets) for key in df.index], names = datasets)
    df = df.iloc[:,0]
    df.name = 'value'

    # Create an UpSet plot
    plot(df, sort_by='cardinality')

    # Save the UpSet plot as a PDF file
    plt.savefig(f"upset_plot_{mode}.pdf", format="pdf")

    # Show the plot (optional)
    plt.show()
    clean_plots()


## 
dfdrugs = plot_upset(mode='drugs')
plot_upset(mode='proteins')