import sympy as sm
import sympy.physics.mechanics as me
from matplotlib.collections import LineCollection
me.init_vprinting()

q1, q2, q3, q4 = me.dynamicsymbols('q1:5')
u1, u2, u3, u4 = me.dynamicsymbols('u1:5')
T = me.dynamicsymbols('T')

l1, l2, l3, l4, l5 = sm.symbols('l1:6')
m2, m3, m5, g = sm.symbols('m2, m3, m5, g')

N = me.ReferenceFrame('N')
A = N.orientnew('A', 'Axis', (q1, N.z))
B = A.orientnew('B', 'Axis', (q2, -A.z))
C = B.orientnew('C', 'Axis', (q3, -B.z))
D = B.orientnew('D', 'Axis', (q4, -B.z))

P1 = me.Point('P1')
P2 = P1.locatenew('P2', l1*A.x)
P3 = P2.locatenew('P3', l2*B.x)
P4 = P3.locatenew('P4', l3*C.x)
P5 = P3.locatenew('P5', l5*D.x)

loop = P4.pos_from(P1) - l4*N.x
config_con1 = loop.dot(N.x).simplify()
config_con2 = loop.dot(N.y).simplify()
qdots = {q.diff(): u for q, u in zip((q1, q2, q3, q4), (u1, u2, u3 ,u4))}

A.set_ang_vel(N, u1*N.z)
B.set_ang_vel(A, -u2*A.z)
C.set_ang_vel(B, -u3*B.z)
D.set_ang_vel(B, -u4*B.z)

P1.set_vel(N, 0)
P1.vel(N)
P2.v2pt_theory(P1, N, A)
P3.v2pt_theory(P2, N, B)
P5.v2pt_theory(P3, N, D)

mot_con = P5.vel(N).dot(D.y).simplify()
t = me.dynamicsymbols._t
config_con1_dot = config_con1.diff(t).subs(qdots)
config_con2_dot = config_con2.diff(t).subs(qdots)
P4.v2pt_theory(P3, N, C).dot(N.x).simplify()
P4.v2pt_theory(P3, N, C).dot(N.y).simplify()
A.ang_acc_in(N)
P2.acc(N)
P3.acc(N)
P5.acc(N)
particle2 = me.Particle('P2', P2, m2)
particle3 = me.Particle('P3', P3, m3)
particle5 = me.Particle('P5', P5, m5)
particles = [particle2, particle3, particle5]
loads = [(P2, -m2*g*N.y),
         (P3, -m3*g*N.y),
         (P5, -m5*g*N.y),
         (A, T*N.z)]

kane = me.KanesMethod(N, # inertial reference frame
                      (q1, q4), # only two independent generalized coordinates
                      (u1,), # only one independent generalized speed
                      kd_eqs=[qd - u for qd, u in qdots.items()], # q' = u for all coordinates
                      q_dependent=(q2, q3), # two depdendent coordinates from the kinematic loop
                      configuration_constraints=[config_con1,  # two kinematic loop config constraints
                                                 config_con2],
                      u_dependent=(u2, u3, u4), # dependent generalized speeds
                      velocity_constraints=[mot_con,  # nonholonomic constraint
                                            config_con1_dot,  # two holonomic motion constraints
                                            config_con2_dot],
                      # acc constraints are required to ensure all qdots are properly substituted (will not be required in SymPy 1.6)
                      acceleration_constraints=[mot_con.diff(t).subs(qdots),
                                                config_con1_dot.diff(t).subs(qdots),
                                                config_con2_dot.diff(t).subs(qdots)])
fr, frstar = kane.kanes_equations(particles, loads=loads)
zero = fr + frstar
zero.shape
print(me.find_dynamicsymbols(zero))
print(kane.q)
print(kane.u)
print(kane.kindiffdict())
print(me.find_dynamicsymbols(kane.mass_matrix_full))
print(me.find_dynamicsymbols(kane.forcing_full))

import numpy as np
from scipy.optimize import fsolve
from pydy.system import System

sys = System(kane)
l1_val = 1.0
l2_val = 1.0
l3_val = 1.0
l4_val = 1.0
l5_val = 3.0

sys.constants = {g: 9.81,
                 l1: l1_val,
                 l2: l2_val,
                 l3: l3_val,
                 l5: l5_val,
                 m2: 1.0,
                 m3: 1.0,
                 m5: 1.0}

q1_0 = np.deg2rad(30.0)
eval_config_con = sm.lambdify((q1, q2, q3, l1, l2, l3, l4),
                              sm.Matrix([config_con1, config_con2]))
eval_config_con_fsolve = lambda x, q1, l1, l2, l3, l4: np.squeeze(eval_config_con(q1, x[0], x[1], l1, l2, l3, l4))
q2_0, q3_0 = fsolve(eval_config_con_fsolve, np.ones(2), args=(q1_0, l1_val, l2_val, l3_val, l4_val))

sys.initial_conditions = {q1: q1_0,
                          q2: q2_0,
                          q3: q3_0,
                          q4: np.deg2rad(-90.0),
                          u1: 0.0,
                          u2: 0.0,
                          u3: 0.0,
                          u4: 0.0}
duration = 4.0
fps = 60.0
sys.times = np.linspace(0.0, duration, num=int(fps*duration))
def step_pulse(x, t):
    if t < 1.0:
        T = -100.0
    else:
        T = 0.0
    return np.array([T])
sys.specifieds = {T: step_pulse}
x = sys.integrate()
import matplotlib.pyplot as plt
fig, axes = plt.subplots(4, 2, sharex=True)
fig.set_size_inches(10, 10)

for i, (xi, ax, s) in enumerate(zip(x.T, axes.T.flatten(), sys.states)):
    ax.plot(sys.times, np.rad2deg(xi))
    title = sm.latex(s, mode='inline')
    ax.set_title(title)
    if 'q' in title:
        ax.set_ylabel('Angle [deg]')
    else:
        ax.set_ylabel('Angular Rate [deg/s]')

axes[3, 0].set_xlabel('Time [s]')
axes[3, 1].set_xlabel('Time [s]')

plt.tight_layout()
config_constraint_vals = eval_config_con(x[:, 0],  # q1
                                         x[:, 2],  # q2
                                         x[:, 3],  # q3,
                                         l1_val, l2_val, l3_val, l4_val).squeeze()

fig, ax = plt.subplots()
ax.plot(sys.times, config_constraint_vals.T)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Configuration Constraint Value [m]')
eval_motion_con = sm.lambdify((q1, q2, q3, q4, u1, u2, u3, u4, l1, l2, l3, l4, l5),
                              sm.Matrix([mot_con, config_con1_dot, config_con2_dot]))
motion_constraint_vals = eval_motion_con(x[:, 0],  # q1
                                         x[:, 2],  # q2
                                         x[:, 3],  # q3
                                         x[:, 1],  # q4
                                         x[:, 4],  # u1
                                         x[:, 5],  # u2
                                         x[:, 6],  # u3
                                         x[:, 7],  # u4
                                         l1_val, l2_val, l3_val, l4_val, l5_val).squeeze()
fig, ax = plt.subplots()
ax.plot(sys.times, motion_constraint_vals.T)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Motion Constraint Value [m/s]')
ax.legend(['Nonholonomic', 'Holonomic', 'Holonomic'])

from matplotlib.animation import FuncAnimation
P2.pos_from(P1).express(N).simplify()
P3.pos_from(P2).express(N).simplify()
P5.pos_from(P3).express(N).simplify()
q1_vals = x[:, 0]
q2_vals = x[:, 2]
q3_vals = x[:, 3]
q4_vals = x[:, 1]

p2_xy = np.array([l1_val*np.cos(q1_vals),
                  l1_val*np.sin(q1_vals)])

p3_xy = p2_xy + np.array([l2_val*np.cos(q1_vals - q2_vals),
                          l2_val*np.sin(q1_vals - q2_vals)])

p5_xy = p3_xy + np.array([l5_val*np.cos(-q1_vals + q2_vals + q4_vals),
                          -l5_val*np.sin(-q1_vals + q2_vals + q4_vals)])

fig, ax = plt.subplots()
fig.set_size_inches((6, 5))
points = np.array([[0.0, 0.0], [p2_xy[:, 0][0], p2_xy[:, 0][1]], [p3_xy[:, 0][0], p3_xy[:, 0][1]], 
                   [p5_xy[:, 0][0], p5_xy[:, 0][1]], [p3_xy[:, 0][0], p3_xy[:, 0][1]], [l4_val, 0.0]])
segments = [points[i:i+2] for i in range(len(points)-1)]
colors = ['red', 'green', 'blue', 'purple', 'orange']

# line, = ax.plot([0.0, p2_xy[0, 0], p3_xy[0, 0], p5_xy[0, 0], p3_xy[0, 0], l4_val],
#                 [0.0, p2_xy[1, 0], p3_xy[1, 0], p5_xy[1, 0], p3_xy[1, 0], 0.0])
lc = LineCollection(segments, colors=colors, linewidths=2)
ax.add_collection(lc)

title = 'Time = {:0.1f} seconds'
ax.set_title(title.format(0.0))
ax.set_ylim((-1.0, 5.0))
ax.set_xlim((-2.0, 3.0))
ax.set_aspect('equal')

def update(i):
    points = np.array([[0.0, 0.0], [p2_xy[:, i][0], p2_xy[:, i][1]], [p3_xy[:, i][0], p3_xy[:, i][1]], 
                       [p5_xy[:, i][0], p5_xy[:, i][1]], [p3_xy[:, i][0], p3_xy[:, i][1]], [l4_val, 0.0]])
    segments = [points[j:j+2] for j in range(len(points)-1)]
    lc.set_segments(segments)  # Update the segments in the LineCollection
    ax.set_title('Time = {:0.1f} seconds'.format(sys.times[i]))  # Update the title
    return lc,

ani = FuncAnimation(fig, update, save_count=len(sys.times)-1)
ani.save('./pydy_test/1.mp4',fps=fps)
plt.show()

pass