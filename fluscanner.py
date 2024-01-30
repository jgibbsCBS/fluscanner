import sys
import os

SeqName=""

ListOfStrains = []
Cal04="MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADTYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI"

HA2Offset=343

def addStrain(SeqName,line,result):
    HA1changes=0
    myStrain=[]
    myStrain.append(SeqName)
    myStrain.append(str(result))
    for x in (3,6,10,57,60,107,109,110,112,113,115,120,138,149,163):
        AA = Cal04[HA2Offset + x]
        AA2 = line[HA2Offset + x]
        if(AA != AA2 and AA2 != "-"):
            myStrain.append(AA2)
        else:
            myStrain.append("")
    myStrain.append("")
    for x in (41,45,61,145,156,224,309):
        AA = Cal04[x-1]
        AA2 = line[x-1]
        if(AA != AA2 and AA2 != "-"):
            myStrain.append(AA2)
        else:
            myStrain.append("")
    myStrain.append("")

    for x in range(1,345):
        AA= line[x-1]
        Cal04AA= Cal04[x-1]
        if(AA != Cal04AA and AA != "-"):
            HA1changes=HA1changes+1
            mutation=Cal04AA + str(x) + AA
            myStrain.append(mutation)
    if (HA1changes > 0):
        myStrain[25]=str(HA1changes)
    ListOfStrains.append(myStrain)

def addHeader1():
    myLine=[]
    myLine.append("Strain")
    myLine.append("stem changes")
    for x in (3,6,10,57,60,107,109,110,112,113,115,120,138,149,163):
        myLine.append(str(x))
    myLine.append("HA1 changes")
    for x in (41,45,61,145,156,224,309):
        myLine.append(str(x))
    myLine.append("Total HA1")
    ListOfStrains.append(myLine)

def addHeader2():
    myLine=[]
    myLine.append("Cal04")
    myLine.append("")
    for x in (3,6,10,57,60,107,109,110,112,113,115,120,138,149,163):
        start=HA2Offset+x
        AA = Cal04[HA2Offset + x]
        myLine.append(AA)
    myLine.append("")
    for x in (41,45,61,145,156,224,309):
        offset = x - 1
        AA = Cal04[offset]
        myLine.append(AA)
    ListOfStrains.append(myLine)


def handleLine(line,SeqName):
    HA2results=0
    totalresults=0
    codon=0

    for x in (3,6,10,57,60,107,109,110,112,113,115,120,138,149,163):
        AA = Cal04[HA2Offset + x]
        AA2 = line[HA2Offset + x]
        if(AA != AA2 and AA2 != "-"):
            HA2results=HA2results+1
    totalresults=HA2results
    for x in (41,45,61,145,156,224,309):
        AA = Cal04[x-1]
        AA = Cal04[x-1]
        AA2 = line[x-1]
        if (AA != AA2 and AA2 != "-"):
            totalresults = totalresults + 1
    if(totalresults > 0 and totalresults < 10):
        addStrain(SeqName,line,HA2results)

        


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
        if line[0] == ">":
            SeqName=line[1:-1]
        else:
            handleLine(line,SeqName)
    
    print ("ListOfStrains:")
    for x in ListOfStrains:
        print(x)

main()


