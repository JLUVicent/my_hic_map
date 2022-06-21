# 读取txt文件
def read_txt(txt):
    with open(txt, 'r') as f:
        lines = f.readlines()
        for line in lines:
            dataSourceList = data_to_list(line)
    return dataSourceList


dataSourceList = []
FinalDataList = []
# 数据存入列表中


def data_to_list(line):
    dataSourceList.append(line.split("\t"))
    return dataSourceList

# 获取最大值对应的索引值


def get_max_index(numlist, list_str):
    for j in range(len(list_str)):
        numlist.append(float(list_str[j][2]))
    maxindex = numlist.index(max(numlist))
    return maxindex


# t = 0
# 处理dataSourceList中的数据
def process_list(dataSourceList):
    intermediateList = []
    numlist = []
    for i in range(len(dataSourceList)):
        if(i == 0):
            intermediateList.append(dataSourceList[i])
        else:
            if(dataSourceList[i][0] == dataSourceList[i-1][0]):
                intermediateList.append(dataSourceList[i])
                if(i == len(dataSourceList)-1):
                    maxindex = get_max_index(numlist, intermediateList)
                    FinalDataList.append(intermediateList[maxindex])
            elif(dataSourceList[i][0] != dataSourceList[i-1][0] or i == len(dataSourceList)-1):
                maxindex = get_max_index(numlist, intermediateList)
                FinalDataList.append(intermediateList[maxindex])
                # 将intermediateList清空,numlist置空
                numlist = []
                intermediateList = []
                intermediateList.append(dataSourceList[i])
    # FinalDataList.append(intermediateList[0])
    return FinalDataList

# 将列表中的元素写入txt文件


def write_list(final_list):
    with open("vicent616_blast_out.txt", "w") as f:
        for i in range(len(final_list)):
            # 写入全部文件
            # f.write(str(my_list[i][0]) + "\t" + str(my_list[i][1]) +
            #         "\t" + str(my_list[i][2]) + "\t" + str(my_list[i][3]))
            # 写入第一列和第四列
            f.write(str(final_list[i][0]) + "\t" + str(final_list[i][3]))


# 程序入口
if __name__ == '__main__':
    # dataSourceList将读入的文件存入一个列表中，
    dataSourceList = read_txt("/media/ubuntu/conda/vicent/HiCzin/合成酵母/megahit_out/vicent616_blast")
    # 加上最后一行，解决源文件最后一行读不到的问题
    add_help = [0, 0, 0, 0]
    dataSourceList.append(add_help)
    # 存储最终处理后的数据
    FinalDataList = process_list(dataSourceList)
    # 写入txt文件
    write_list(FinalDataList)
    print("处理完成")
