
def process_line(raw_line):
    '''
    >>> process_line('scaffold123 159 snv1 G C 30 PASS NS=7;DP=38;AF=0.071;AA=G GT:DP 0/0:6 0/0:3 0/0:6 0/0:4 0/1:5 0/0:9 0/0:5')
    'G 5 8 C 1 0 scaffold123 159'
    
    >>> process_line('scaffold123 15934 snv20 A C 30 PASS NS=7;DP=56;AF=0.214;AA=A GT:DP 0/0:5 0/0:10 0/0:9 0/1:12 0/0:8 0/1:6 0/1:6')
    'A 6 5 C 0 3 scaffold123 15934'
    '''

    raw_line = raw_line.split()
    allele1 = raw_line[3]
    allele2 = raw_line[4]
    chrom = raw_line[0]
    position = raw_line[1]
    pops_allele1 = [0, 0]
    pops_allele2 = [0, 0]

    def getTan(pair):
        index1 = int(pair[0])
        index2 = int(pair[2])

        if index1 == 0:
            pops_allele1[0] += 1
        else:
            pops_allele2[0] += 1

        if index2 == 0:
            pops_allele1[0] += 1
        else:
            pops_allele2[0] += 1

    def getNam(pair):
        index1 = int(pair[0])
        index2 = int(pair[2])

        if index1 == 0:
            pops_allele1[1] += 1
        else:
            pops_allele2[1] += 1

        if index2 == 0:
            pops_allele1[1] += 1
        else:
            pops_allele2[1] += 1

    Tans = [10, 11, 13]
    Nams = [9, 12, 14, 15]

    for index in Tans:
        getTan(raw_line[index])

    for index in Nams:
        getNam(raw_line[index])

    return allele1 + ' ' + str(pops_allele1[0]) + ' ' + str(pops_allele1[1]) + ' ' + allele2 + ' ' +\
           str(pops_allele2[0]) + ' ' + str(pops_allele2[1]) + ' ' + chrom + ' ' + position
# with open('aju-ANGSD-all_updated.vcf') as input, open('aju-ANGSD-dadi.txt', 'w') as output:
#     for i in range(10):
#         line = input.readline()
#     output.write('#Allele1 Tan Nam Allele2 Tan Nam Chromosome Position\n')
#
#     for line in input:
#         line = process_line(line)
#         output.write(line + '\n')


def isTriplet(str):
    return str.isalpha()

with open('result_fin.vcf') as input, open('fin_result_dadi.txt', 'w') as output:
    # for i in range(10):
    #     line = input.readline()
    output.write('Cheetah Catan Allele1 Tan Nam Allele2 Tan Nam Chromosome Position\n')

    for line in input:
        splited = line.split()
        cheetah, catan = splited[-2:]
        print(cheetah, catan)
        if len(cheetah) != 3 or len(catan) != 3 :
            continue
        if not isTriplet(cheetah) or not isTriplet(catan):
            continue
        line = process_line(line)
        output.write(cheetah + ' ' + catan + ' ' + line + '\n')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
