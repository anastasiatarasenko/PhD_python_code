# coding: utf-8

import datetime

__author__ = 'Alexander Myasoedov'
__email__ = 'mag@rshu.ru'
__created__ = datetime.datetime(2016, 5, 11)

#
# Слушатель вызывает Processor / BaseHandler для чтения и т.д. с заданными параметрами:
# Возможность изменять:
#     фильтр по имени файла
#     входной/выходной пути
#     временной диапазон
#     разрешение
#     район по широтам/долготам
#     максимальный размер снимка


class Listener:
    def __init__(self):
        pass
