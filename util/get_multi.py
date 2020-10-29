import pandas as pd 

utr5_pool1 = pd.read_excel(open('UTR5.xlsx','rb'),sheet_name='Pool1',index_col=0) 
utr5_pool2 = pd.read_excel(open('UTR5.xlsx','rb'),sheet_name='Pool2',index_col=0)


utr5_pool1_3 = utr5_pool1 
for i in range(2):
    utr5_pool1_3 = utr5_pool1_3.append(utr5_pool1,ignore_index=True) 

utr5_pool2_44 = utr5_pool2
for i in range(43):
    utr5_pool2_44 = utr5_pool2_44.append(utr5_pool2,ignore_index=True)

writer = pd.ExcelWriter('utr5_multiple.xlsx', engine='xlsxwriter')
utr5_pool1_3.to_excel(writer,sheet_name='Pool1')
utr5_pool2_44.to_excel(writer,sheet_name='Pool2') 
writer.save()

