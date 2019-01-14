import requests
from lxml import etree
import os
import time

def getHtml(url):
    '''
    get the corresponding response of the url and decode it
    :param url: url
    :return: decoded response
    '''
    html = requests.get(url).content
    selector = etree.HTML(html)
    return selector

def getContent(htm, xpathStr):
    '''
    get the content corresponding the partten of xpathStr(title paper id etc.)
    :param htm: response
    :param xpathStr: xpath partten
    :return: content
    '''
    selector = htm
    content = selector.xpath(xpathStr)
    return content

def getDownPdf(pdf_id,pdf_title,url0):
    '''
    download paper choosen
    :param pdf_id: paper's id
    :param pdf_title: paper's title
    :param url0: url where we download the paper
    :return: None
    '''

    print('downloading '+pdf_title+' ....')
    os.chdir('/Users/shiwuzhang/ASTRO/paper')
    print(pdf_id)
    url=url0+pdf_id+'.pdf'
    print(url)
    header={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.110 Safari/537.36'}
    paper=open(pdf_title+'.pdf','wb')
    response=requests.get(url,headers=header)
    paper.write(response.content)
    paper.close()
    print('Done!')


def PaperCheckDonwload(url,xp1,xp2):
    '''
    check recent paper and choose what kind of paper to download
    :param url: address
    :param xp1: partten to select paper's id
    :param xp2: partten to select paper's title
    :return: None
    '''


    htm0 = getHtml(url)
    cons_id = getContent(htm0, xp1)  # get pdfs' href
    cons_title=[i[:-1] for i in getContent(htm0, xp2) if i!='\n']# get papers' title

    print(cons_id)
    print(cons_title)
    selection=input('Do you want to download all?')
    if selection=='yes':
        for i in range(len(cons_title)):
            getDownPdf(cons_id[i],cons_title[i],'https://www.arxiv.org')
    else:
        a=input('what kind of paper do you want?(only for one words)')
        paper=[(cons_id[i],cons_title[i]) for i in range(len(cons_id)) if a in cons_title[i]]
        print(paper)
        for p in paper:
            getDownPdf(p[0],p[1],'https://www.arxiv.org')

if __name__=='__main__':
    url='https://arxiv.org/list/astro-ph/pastweek?show=361'#update url before get papers
    xp1='//dt//*[@class="list-identifier"]//a[2]//@href'
    xp2='//div[@class="list-title mathjax"]/text()'
    PaperCheckDonwload(url,xp1,xp2)
