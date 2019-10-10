import sys

SYMBOLS = {'}':'{', ']':'[', ')':'('}
def match_brackets(s):
    bracket = []
    for i in s:
        if i in SYMBOLS.values():
             bracket.append(i)
        elif i in SYMBOLS.keys():
             if bracket and bracket[-1] == SYMBOLS[i]:
                 bracket.pop()
             else:
                 return False
    return not bracket
        

print(match_brackets("{{([({})]}}"))
