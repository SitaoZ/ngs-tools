#coding:utf-8
import re,sys
import xlrd



data = xlrd.open_workbook('huhaiyan.xlsx')
out = open('haiyan.xlsx','w')

sheet = data.sheet_by_name('Sheet1') #行数
print (sheet.nrows)
#print (sheet.name,sheet.ncols,sheet.nrows)
#for rownum in range(1,sheet.nrows):
#   out.writelines("\t".join(sheet.row_values(rownum))+"\n")
#out.close()
