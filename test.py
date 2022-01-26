import shelve
import matplotlib.pyplot as plt
# with shelve.open("/home/DUAN/methylation/data_files/ucsc") as db:
#     dic = db["ucsc"]
# # new_dic = {i:dic[i].values() for i in dic}
# # print(new_dic)

# new_list = list(dic.keys())[:9]+list(dic.keys())[10:25]
# print(new_list)
# # print(list(dic["chr1"].values()))
# new_dic = {i:list(dic[i].values()) for i in new_list}

# with shelve.open("/home/DUAN/methylation/data_files/dic") as dbs:
#     dbs["dic"] = new_dic

# with shelve.open("/home/DUAN/methylation/data_files/dic") as dbs:
#     new_dic = dbs["dic"]
# print(new_dic.keys())
# # def func(listTemp, n):
# #     for i in range(0, len(listTemp), n):
# #         yield listTemp[i:i + n]
# # a=len(range(247754845, 247760556))
# # b = 247754845+620
# # print(a)
# # print(list(range(247754845, 247760556))[0])
# # print((b-list(range(247754845, 247760556))[0])/a)

# help_dict = {"-2kb":0,"10%":0,"20%":0,"30%":0,"40%":0,"50%":0,"60%":0,"70%":0,"80%":0,"90%":0,"100%":0,"+2kb":0}

# plt.plot(x,y)
# plt.show()
# plt.savefig("/home/DUAN/plt.png")
dic = {'-2kb': 0, '10%': 5934, '20%': 5749, '30%': 5679, '40%': 5613, '50%': 5578, '60%': 5666, '70%': 5972, '80%': 5621, '90%': 5796, '100%': 6058, '+2kb': 0}

y = list(dic.values())[1:11]
x = list(dic.keys())[1:11]
plt.plot(x,y)
plt.xlabel('gene position sliced to ten equal parts')
plt.ylabel('methylation count in certain area')
plt.show()
plt.savefig("/home/DUAN/plt.png")