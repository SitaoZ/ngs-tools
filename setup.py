# -*- coding: UTF-8 -*-
import os

import setuptools

setuptools.setup(
    name='ngs-tools',
    version='2019.03.31',
    keywords='demo',
    description='A demo for python packaging.',
    long_description=open(
        os.path.join(
            os.path.dirname(__file__),
            'README.rst'
        )
    ).read(),
    author='zhusitao1990',  # 替换为你的Pypi官网账户名
    author_email='zhusitao1990@163.com',  # 替换为你Pypi账户名绑定的邮箱

    url='https://github.com/zhusitao1990/ngs-tools',  # 这个地方为github项目地址，貌似非必须
    packages=setuptools.find_packages(),
    license='MIT'
)