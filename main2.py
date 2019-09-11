"""
Version: 1.1

The Script retruns all possible cuts and how likely they are for a specific protease and sequence.



"""


import MySQLdb
import data
import operator

CODES = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
         "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
         "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U",
         "Pyl": "O", "Xaa": "X", "Asx": "B", "Glx": "Z", "Xle": "J", "-": "-"}

def _write_data(data): #writes data into data.py file
    f = open("data.py", "w")
    f.write("data = " + data)
    f.close()

def translate(_3letter): #translates 3-letter-code into 1-letter-code
    if _3letter in CODES.keys():
        return CODES[_3letter]
    else:
        return "X"

def print_results(r): #sorts results according to score, prints them
    sorted_r = sorted(r.items(), key=operator.itemgetter(1))
    print(sorted_r[::-1])

def score(data_list): #assigns each cut in data_list a score how likely that cut is
    scoredict = {}
    for k in data_list.keys():
        scoredict[k] = data_list[k][3] + data_list[k][4]
    return scoredict

def get_data(protease, data): #checks if data is in data.py, adds data to data.py if it's not, returns data
    if protease in data.data:
        print("FOUND IN DATA")
    else:
         table = calc_table(protease)
         data.data[protease] = table
         _write_data(str(data.data))
         print("ADDED TO DATA")
    return data.data[protease]

def calc_table(protease): #gets data from sql, returns a dictionary of how many times the protease cuts after each specific amino-acid
    conn = MySQLdb.connect(host='localhost',user='root',passwd='SH01m7n9i_.', db='merops121')
    cursor = conn.cursor()
    command = "SELECT Site_P4, Site_P3, Site_P2, Site_P1, Site_P1prime, Site_P2prime, Site_P3prime, Site_P4prime FROM Substrate_search WHERE Protease = '{0}'".format(protease)
    cursor.execute(command)
    data_list = cursor.fetchall()
    if data_list: #checks if sql returned something
        p_data = {} #writes data from sql to a dictionary
        for v in CODES.values():
            p_data[v] = [0, 0, 0, 0, 0, 0 , 0, 0]
        for l in data_list:
            for i in range(8):
                p_data[translate(l[i])][i] += 1
        #print(data_list)
    else:
        raise ValueError("Protease named {0} not found".format(protease))
    return p_data

def lookup(protease, sequence): #main function 
    sequence = "---" + sequence + "---"
    protease_data = get_data(protease, data) #gets data
         
    cut_dict = {} #creates a dict of all possible cuts and a list of how many times the protease cuts the sequence after that amino-acid
    for i in range(3, len(sequence)-4):
        p_string = sequence[i-3:i+5]
        l = []
        for j in range(len(p_string)):
            l.append(protease_data[p_string[j]][j])
        cut_dict[p_string] = l
    print(cut_dict)
        
    scoredict = score(cut_dict) #scores each cut according to likelihood
    #print_results(scoredict) #prints the rusults
    return scoredict

if __name__ == "__main__":
    #demo
    data_list = lookup('kallikrein-related peptidase 3', "MYREWVVVNVFMMLYVQLVQGSSNEHGPVKRSSQSTLERSEQQIRAASSLEELLRITHSEDWKLWRCRLRLKSFTSMDSRSASHRSTRFAATFYDIETLKVIDEEWQRTQCSPRETCVEVASELGKSTNTFFKPPCVNVFRCGGCCNEESLICMNTSTSYISKQLFEISVPLTSVPELVPVKVANHTGCKCLPTAPRHPYSIIRRSIQIPEEDRCSHSKKLCPIDMLWDSNKCKCVLQEENPLAGTEDHSHLQEPALCGPHMMFDEDRCECVCKTPCPKDLIQHPKNCSCFECKESLETCCQKHKLFHPDTCSCEDRCPFHTRPCASGKTACAKHCRFPKEKRAAQGPHSRKNP")
    
