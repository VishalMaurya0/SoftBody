🧩 Soft Body Physics Engine (C++ / SFML)








A real-time 2D Soft Body Physics Simulation built from scratch using C++ and SFML.

This project implements:

Spring-based deformation

Torque & angular physics

Impulse-based collision resolution

Soft body vs soft body interaction

Custom edge projection collision detection

🎥 Demo


Download this: [Soft Body Demo](media/demo.mp4)

Short version: [Soft Body Demo](media/clip.mp4)

🚀 Features

✅ Spring-mass soft body system (Hooke’s Law)
✅ Center of Mass & Moment of Inertia calculation
✅ Angular velocity & torque simulation
✅ Soft body ↔ soft body collision
✅ Boundary collision with restitution
✅ Friction simulation
✅ Real-time FPS counter
✅ Keyboard interaction controls

💥 Collision Handling

Soft Body vs Boundary

Velocity reflection

Energy loss (restitution)

Torque from impulse

Soft Body vs Soft Body

Closest edge detection

Projection onto edge

Penetration correction

Impulse resolution

Friction handling

🎮 Controls
Key	Action
⬅️ Left Arrow	Apply left force
➡️ Right Arrow	Apply right force
⬆️ Up Arrow	Apply negative torque
⬇️ Down Arrow	Apply positive torque
Spacebar	Apply upward impulse


Core Components

Point → Mass particle

Spring → Constraint between two points

Edge → Used for collision detection

Shape → Soft or rigid body container

CalculatePhysics() → Main physics update

HandleSoftBodyCollision() → Soft body interaction

🛠️ Installation & Build
1️⃣ Install SFML

Download from:
https://www.sfml-dev.org/download.php

Link against:

sfml-graphics

sfml-window

sfml-system

2️⃣ Compile (g++ example)
g++ main.cpp -o softbody -lsfml-graphics -lsfml-window -lsfml-system


👨‍💻 Author

Vishal Maurya

If you found this interesting, consider ⭐ starring the repository.
