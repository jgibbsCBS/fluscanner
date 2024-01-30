import sys
import os

SeqName=""

#ListOfStrains = [["dummyName1",1,24,"NSS"],["dummyName2",1,25,"NTT"]]
#ListOfStrains = [["A/Antarctica/Dummy1/2022","A/Antarctica/Dummy2/2022"]]
ListOfStrains = []
#Cal04="---------------------------------ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCGAACAATTCAACAGACACTGTAGACACAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCTTCTAGAAGACAAGCATAACGGGAAACTATGCAAACTAAGAGGGGTAGCCCCATTGCATTTGGGTAAATGTAACATTGCTGGCTGGATCCTGGGAAATCCAGAGTGTGAATCACTCTCCACAGCAAGCTCATGGTCCTACATTGTGGAAACACCTAGTTCAGACAATGGAACGTGTTACCCAGGAGATTTCATCGATTATGAGGAGCTAAGAGAGCAATTGAGCTCAGTGTCATCATTTGAAAGGTTTGAGATATTCCCCAAGACAAGTTCATGGCCCAATCATGACTCGAACAAAGGTGTAACGGCAGCATGTCCTCATGCTGGAGCAAAAAGCTTCTACAAAAATTTAATATGGCTAGTTAAAAAAGGAAATTCATACCCAAAGCTCAGCAAATCCTACATTAATGATAAAGGGAAAGAAGTCCTCGTGCTATGGGGCATTCACCATCCATCTACTAGTGCTGACCAACAAAGTCTCTATCAGAATGCAGATACATATGTTTTTGTGGGGTCATCAAGATACAGCAAGAAGTTCAAGCCGGAAATAGCAATAAGACCCAAAGTGAGGGATCAAGAAGGGAGAATGAACTATTACTGGACACTAGTAGAGCCGGGAGACAAAATAACATTCGAAGCAACTGGAAATCTAGTGGTACCGAGATATGCATTCGCAATGGAAAGAAATGCTGGATCTGGTATTATCATTTCAGATACACCAGTCCACGATTGCAATACAACTTGTCAAACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCGATCACAATTGGAAAATGTCCAAAATATGTAAAAAGCACAAAATTGAGACTGGCCACAGGATTGAGGAATATCCCGTCTATTCAATCTAGAGGCCTATTTGGGGCCATTGCCGGTTTCATTGAAGGGGGGTGGACAGGGATGGTAGATGGATGGTACGGTTATCACCATCAAAATGAGCAGGGGTCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACGAGATTACTAACAAAGTAAATTCTGTTATTGAAAAGATGAATACACAGTTCACAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGTTGGTTCTATTGGAAAATGAAAGAACTTTGGACTACCACGATTCAAATGTGAAGAACTTATATGAAAAGGTAAGAAGCCAGCTAAAAAACAATGCCAAGGAAATTGGAAACGGCTGCTTTGAATTTTACCACAAATGCGATAACACGTGCATGGAAAGTGTCAAAAATGGGACTTATGACTACCCAAAATACTCAGAGGAAGCAAAATTAAACAGAGAAGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCAGTTTCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATTTAA"
Cal04="MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"

#HA2Offset=1064
HA2Offset=343

def addStrain(SeqName,line,result):
    HA1changes=0
    myStrain=[]
    myStrain.append(SeqName)
    myStrain.append(str(result))
    for x in (3,6,10,45,57,60,107,109,112,113,115,120,149,163):
        AA = Cal04[HA2Offset + x]
        AA2 = line[HA2Offset + x]
        if(AA != AA2 and AA2 != "-"):
            myStrain.append(AA2)
        else:
            myStrain.append("")
    myStrain.append("")
    for x in range(1,345):
        AA= line[x-1]
        Cal04AA= Cal04[x-1]
        # if(Cal04AA == 'X'):
        #     print("position " + str(offset))
#        print (line[offset:offset+3])
#        print (Cal04[offset:offset+3])
#        print(str(x) + " " + str(offset))

        if(AA != Cal04AA and AA != "-"):
            HA1changes=HA1changes+1
            mutation=Cal04AA + str(x) + AA
            myStrain.append(mutation)
    myStrain[16]=str(HA1changes)
    ListOfStrains.append(myStrain)

def addHeader1():
    myLine=[]
    myLine.append("Strain")
    myLine.append("stem changes")
    for x in (3,6,10,45,57,60,107,109,112,113,115,120,149,163):
        myLine.append(str(x))
    myLine.append("HA1 changes")
    ListOfStrains.append(myLine)

def addHeader2():
    myLine=[]
    myLine.append("Cal04")
    myLine.append("")
    for x in (3,6,10,45,57,60,107,109,112,113,115,120,149,163):
        offset=x-1
        #AA= aaFromNT(Cal04[HA2Offset+offset:HA2Offset+offset+3])
        start=HA2Offset+x
        marker=Cal04[start:6]
        AA = Cal04[HA2Offset + x]
        AA2 = Cal04[HA2Offset + x+1]
        AA3 = Cal04[HA2Offset + x+2]
        AA4 = Cal04[HA2Offset + x+3]
        AA5 = Cal04[HA2Offset + offset+4]
        myLine.append(AA)
    ListOfStrains.append(myLine)


def handleLine(line,SeqName):
    result=0
    codon=0
#     while(codon<566):
#         start=1
#         AA = line[start+codon]
# #        print (line[start+offset:start+offset+15])
#         codon=codon+1
#         if AA == "X":
# #            print("dumping " + SeqName)
#            return

    for x in (3,6,10,45,57,60,107,109,112,113,115,120,149,163):
        #AA= aaFromNT(line[HA2Offset+offset:HA2Offset+offset+3])
        AA = Cal04[HA2Offset + x]
# AA= aaFromNT(Cal04[HA2Offset+offset:HA2Offset+offset+3])
        AA = Cal04[HA2Offset + x]
        AA2 = line[HA2Offset + x]
#        print(x)
#        print(offset)
#        print(AA)
#        print(line[HA2Offset+offset:HA2Offset+offset+15])
#         if aaFromNT(line[HA2Offset+offset:HA2Offset+offset+3]) == 'X':
#             return
        if(AA != AA2 and AA2 != "-"):
            result=result+1
    if(result > 0 and result < 8):
        addStrain(SeqName,line,result)
#        print("appending " + SeqName)

        
def Cal04Residue(codon):
    offset=codon*3-2
    Cal04AA= aaFromNT(Cal04[HA2Offset+offset:HA2Offset+offset+3])
#    print(str(codon) + " " + str(offset))
#    print(Cal04[HA2Offset+offset:HA2Offset+offset+15])
    return(Cal04AA)


def matches(position,AA):
    offset=position*3-2
    Cal04AA= aaFromNT(Cal04[HA2Offset+offset:HA2Offset+offset+3])
    if(Cal04AA == AA):
        return 1
    else:
        return 0


def aaFromNT(codon):
#    print("codon is " + codon)
    if codon == "TTT" or codon == "TTC":
        return "F"
    elif codon == "TTA" or codon == "TTG" or codon[0:2] == "CT":
        return "L"
    elif codon[0:2] == "TC" or codon == "AGT" or codon == "AGC":
        return "S"
    elif codon == "TAT" or codon == "TAC":
        return "Y"
    elif codon == "TGT" or codon == "TGC":
        return "C"
    elif codon == "TGG":
        return "W"
    elif codon[0:2] == "CC":
        return "P"
    elif codon == "CAT" or codon == "CAC":
        return "H"
    elif codon == "CAA" or codon == "CAG":
        return "Q"
    elif codon[0:2] == "CG" or codon == "AGA" or codon == "AGG":
        return "R"
    elif codon == "ATT" or codon == "ATC" or codon == "ATA":
        return "I"
    elif codon == "ATG":
        return "M"
    elif codon[0:2] == "AC":
        return "T"
    elif codon == "AAT" or codon == "AAC":
        return "N"
    elif codon == "AAA" or codon == "AAG":
        return "K"
    elif codon[0:2] == "GT":
        return "V"
    elif codon[0:2] == "GC":
        return "A"
    elif codon == "GAT" or codon == "GAC":
        return "D"
    elif codon == "GAA" or codon == "GAG":
        return "E"
    elif codon[0:2] == "GG":
        return "G"
    else:
#        print("codon is " + codon)
        return "X"


def main():

    print("\nBeginning fluscanner")
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the input file.")
    infile = args[0]
    print("Reading input from %s" % infile)

    if not os.path.isfile(infile):
        raise IOError("Cannot find infile of %s in the current directory." % infile)

    file = open(infile, "r")

    addHeader1()
    addHeader2()

    for line in file:
#        if not line.isspace():
        if line[0] == ">":
            SeqName=line[1:-1]
#            print SeqName
        else:
            handleLine(line,SeqName)
    
    print ("ListOfStrains:")
    for x in ListOfStrains:
        print(x)
#    print(ListOfStrains)

main()


