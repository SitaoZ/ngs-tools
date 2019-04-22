#coding:utf-8
import sys
import xlrd
data = xlrd.open_workbook('AD1_HAU_v1.0_genes2Go.xlsx')
names = data.sheet_names()
out = open('cotton.GO','w')
for name in names:
    sheet = data.sheet_by_name(name) #行数
    #print (sheet.name,sheet.ncols,sheet.nrows)
    for rownum in range(3,sheet.nrows):
       out.writelines("\t".join(sheet.row_values(rownum))+"\n")
out.close()






