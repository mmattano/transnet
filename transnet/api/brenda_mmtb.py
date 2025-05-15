import urllib.request
import pandas as pd
import traceback

# Parameters
Dir_input = "./Resource/"
Dir_output = "./Output/"
FileName_input = "List-cpd (230306).txt"
FileName_output = "KEGG2name_expanded (lower).txt"
url_header = "http://mmtb.tu-bs.de/search/"

# Load Files
compounds_url = "https://rest.kegg.jp/list/compound"

compounds = pd.read_table(compounds_url, sep="\t", header=None)
cols = compounds.columns
compounds.columns = ["compound_id", "name"]
# Expand Names
#
#

new_names = []
for _, row in compounds.iterrows():
    temp_KEGG = row["compound_id"]
    temp_Name = row["name"]
    temp_searchname = temp_Name.replace(" ", "%20")
    temp_url = url_header + temp_searchname
    temp_NewName_set = set()
    temp_NewName_set.add(temp_Name.lower())
    try:
        urllib.request.urlretrieve(temp_url, "./temp.txt")
        f = open("./temp.txt", encoding="utf8")
        text = f.read()
        f.close()
        if "\n    <li>" not in text:
            print(temp_Name + " cannot be found!")
            f_new = open(Dir_output + FileName_output, "a+", encoding="utf8")
            f_new.write(temp_KEGG + "\t" + temp_Name.lower() + "\n")
            f_new.close()
            continue
        else:
            temp_list_1 = text.split("\n    <li>")
            for temp_1 in temp_list_1[1 : len(temp_list_1)]:
                temp_name_new = temp_1.split("</li>")[0]
                temp_name_new = temp_name_new.replace("_", " ")
                temp_name_new = temp_name_new.replace("&#39;", "'")
                temp_name_new = temp_name_new.lower()
                temp_NewName_set.add(temp_name_new)
            temp_NewName_list = list(temp_NewName_set)
            temp_NewName_list.sort()
            new_entry = {
                "KEGG": temp_KEGG,
                "Name": temp_NewName_list,
            }
            new_names.append(new_entry)
            # f_new = open(Dir_output + FileName_output, 'a+', encoding = 'utf8')
            # f_new.write(temp_KEGG + '\t' + ';'.join(temp_NewName_list) + '\n')
            # f_new.close()
            print(temp_KEGG + " is done!")
    except:
        print(temp_Name + " cannot be found!")
        # f_error = open('./error.txt', 'a+', encoding = 'utf8')
        # f_error.write(row)
        # f_error.close()
