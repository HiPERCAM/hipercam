#!/usr/bin/env python

import matplotlib.pyplot as plt

noclear = (
    ('C', True),
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

plt.plot([8.5,16.5,16.5,8.5,8.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([20.5,17.5,17.5,20.5],
         [-0.06,-0.06,0.08,0.08],'--b')

plt.plot([0.5,7.5,7.5,0.5,0.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--r')

plt.xlim(-0.5,len(noclear)-0.3)
plt.ylim(-0.45,0.45)
plt.axis('off')
plt.savefig('noclear.png')


clear = (
    ('C', True),
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

plt.plot([12.5,21.5,21.5,12.5,12.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([26.5,24.5,24.5,26.5],
         [-0.06,-0.06,0.08,0.08],'--b')

plt.plot([0.5,9.5,9.5,0.5,0.5],
         [-0.06,-0.06,0.08,0.08,-0.06],'--r')

plt.xlim(-0.5,len(clear)-0.3)
plt.ylim(-0.45,0.45)
plt.axis('off')
plt.savefig('clear.png')

drift = (
    ('C', True),
    ('E', True),
    ('LS', False),
    ('LD', False),
    ('R', False),
    ('E', False),
    ('LS', False),
    ('LD', False),
    ('R', False),
    ('E', False),
    ('LS', True),
    ('LD', True),
    ('R', True),
    ('E', False),
    ('LS', True),
    ('LD', True),
    ('R', True),
    ('E', False),
    ('LS', True),
    )

fig = plt.figure(figsize=(8,2))
ax = plt.axes([0,0,1,1])

nt = 1
nd = -2 # DRIFT NWINS = 3
for n, (letter, flag) in enumerate(drift):
    x = 1.5*n
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
        plt.text(x,+0.4,'F{:d}'.format(nd),color='b',
                 ha='center',va='center')
        nd += 1

plt.plot([0.8,2.,2.,0.8,0.8],
         [-0.06,-0.06,0.08,0.08,-0.06],'--r')

plt.plot([3.8,8.,8.,3.8,3.8],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([9.8,14.,14.,9.8,9.8],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([15.8,20.,20.,15.8,15.8],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.plot([21.8,26.,26.,21.8,21.8],
         [-0.06,-0.06,0.08,0.08,-0.06],'--b')

plt.xlim(-0.5,len(clear)-0.3)
plt.ylim(-0.45,0.45)
plt.axis('off')
plt.savefig('drift.png')
