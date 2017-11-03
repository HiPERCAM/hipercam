#!/usr/bin/env python

import matplotlib.pyplot as plt

noclear = (
    ('C', True),
    ('W', True),
    ('E', True),
    ('F', False),
    ('R', False),
    ('E', False),
    ('F', False),
    ('R', False),
    ('E', False),
    ('F', True),
    ('R', True),
    ('E', True),
    ('F', False),
    ('R', False),
    ('E', False),
    ('F', False),
    ('R', False),
    ('E', False),
    ('F', True),
    ('R', True),
    ('E', True),
    ('F', False),
    )


fig = plt.figure(figsize=(8,2))
ax = plt.axes([0,0,1,1])

nt = 1
nd = 1
for n, (letter, flag) in enumerate(noclear):
    x = n
    if flag:
        plt.text(x,0,letter,color='g',ha='center',va='center',size=18)
    else:
        plt.text(x,0,letter,color='r',ha='center',va='center', size=18)

    if letter == 'E':
        x -= 0.4
        plt.plot([x,x],[-0.08,-0.3],'k')
        plt.text(x,-0.4,'TS{:d}'.format(nt),ha='center',va='center')
        nt += 1
    elif letter == 'R':
        x += 0.4
        plt.plot([x,x],[0.11,0.3],'k')
        if nd % 3 == 0:
            plt.text(x,+0.4,'F{:d}'.format(nd),color='b',
                     ha='center',va='center')
        else:
            plt.text(x,+0.4,'F{:d}'.format(nd),color='k',
                     ha='center',va='center')
        nd += 1

plt.plot([9.5,17.5,17.5,9.5,9.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([21.5,18.5,18.5,21.5],
         [-0.06,-0.06,0.08,0.08],'--b')

plt.plot([1.5,8.5,8.5,1.5,1.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--r')

plt.xlim(-0.5,len(noclear)-0.3)
plt.ylim(-0.45,0.45)
plt.axis('off')
plt.savefig('noclear.pdf')


clear = (
    ('C', True),
    ('W', True),
    ('E', True),
    ('F', False),
    ('R', False),
    ('W', False),
    ('E', False),
    ('F', False),
    ('R', False),
    ('W', False),
    ('E', False),
    ('F', True),
    ('R', True),
    ('W', True),
    ('E', True),
    ('F', False),
    ('R', False),
    ('W', False),
    ('E', False),
    ('F', False),
    ('R', False),
    ('W', False),
    ('E', False),
    ('F', True),
    ('R', True),
    ('W', True),
    ('E', True),
    ('F', False),
    )

fig = plt.figure(figsize=(8,2))
ax = plt.axes([0,0,1,1])

nt = 1
nd = 1
for n, (letter, flag) in enumerate(clear):
    x = n
    if flag:
        plt.text(x,0,letter,color='g',ha='center',va='center',size=18)
    else:
        plt.text(x,0,letter,color='r',ha='center',va='center', size=18)

    if letter == 'E':
        x -= 0.5
        plt.plot([x,x],[-0.08,-0.3],'k')
        plt.text(x,-0.4,'TS{:d}'.format(nt),ha='center',va='center')
        nt += 1
    elif letter == 'R':
        x += 0.5
        plt.plot([x,x],[0.11,0.3],'k')
        if nd % 3 == 0:
            plt.text(x,+0.4,'F{:d}'.format(nd),color='b',
                     ha='center',va='center')
        else:
            plt.text(x,+0.4,'F{:d}'.format(nd),color='k',
                     ha='center',va='center')
        nd += 1

plt.plot([13.5,22.5,22.5,13.5,13.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([27.5,25.5,25.5,27.5],
         [-0.06,-0.06,0.08,0.08],'--b')

plt.plot([1.5,10.5,10.5,1.5,1.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--r')

plt.xlim(-0.5,len(clear)-0.3)
plt.ylim(-0.45,0.45)
plt.axis('off')
plt.savefig('clear.pdf')
