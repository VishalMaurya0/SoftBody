// FirstC++.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <SFML/Graphics.hpp>
#include <optional> // needed for std::optional
#include <cmath> // needed for sqrt, atan2, sin, cos
#include <chrono>
#include <thread>
#include "SoftBody.h"
#include <iostream>

using sf::Vector2f;
using sf::Color;
using sf::RenderWindow;
using sf::RectangleShape;
using sf::CircleShape;
using sf::VideoMode;
using sf::Event;
using sf::Vertex;
using sf::PrimitiveType;

using std::vector;

struct Point
{
    Vector2f pos;
    Vector2f velocity;
    Vector2f acceleration;
    Vector2f force;
    Vector2f addForce;

    float mass = 1;

    Color color = Color::White;
};

struct Edge
{
    Point* point1;
    Point* point2;

    Color color = Color::White;
};

struct Spring
{
    Point& point1;
    Point& point2;
    float stiffnes;
    float restLength;

    Spring(Point& point1, Point& point2, float stiffnes) : point1(point1), point2(point2), stiffnes(stiffnes)
    {
        restLength = std::sqrt((point2.pos.x - point1.pos.x) * (point2.pos.x - point1.pos.x) + (point2.pos.y - point1.pos.y) * (point2.pos.y - point1.pos.y));
    }
};



class Shape {
    public: vector<Point> points;
    public: bool isSoftBody = true;
    public: float stiffness;
    public: bool staticBody = false;
    public: Color color = Color::White;
    
// Stored Data
    public: vector<Spring> springs;
    public: vector<Edge> edges;

// Dynamic Data
    public: Vector2f centerOfMass;
    public: float angularVelocity = 0.f;
    public: float torque = 0.f;
    public: float momentOfInertia = 0.f;
public: Vector2f addForce = {0.f, 0.f};



    Shape(const vector<Point>& points, const bool isSoftBody, const float stiiffness = 2) {
        this->isSoftBody = isSoftBody;
        this->points = points;
        this->stiffness = stiiffness;
        CalculateCenterOfMass();
        CreateSprings();
        CreateEdges();
    }

public:
    void CalculateCenterOfMass() 
    {
        if (points.empty()) {
            centerOfMass = Vector2f(0.f, 0.f);
            return;
        }

        Vector2f sum(0.f, 0.f);
        for (const auto& point : points) {
            sum += point.pos;
        }
        centerOfMass = sum / static_cast<float>(points.size());


        // Approximate moment of inertia
        momentOfInertia = 0.f;
        for (auto& point : points) {
            Vector2f r = point.pos - centerOfMass;
            float dist2 = r.x * r.x + r.y * r.y;
            momentOfInertia += point.mass * dist2;
        }
    }

public:
    float CalculateMomentOfInertia()
    {
        float momentOfInertia = 0.f;
        for (auto& point : points) {
            Vector2f r = point.pos - centerOfMass;
            float dist2 = r.x * r.x + r.y * r.y;
            momentOfInertia += point.mass * dist2;
        }
        return momentOfInertia;
    }

private:
    void CreateSprings()
    {
        for (size_t i = 0; i < points.size(); i++)
        {
            int pointIndex = i;
            int nextPointIndex = (i + 1) % points.size();
            int nextNextPointIndex = (i + 2) % points.size();

            springs.push_back(Spring(points[pointIndex], points[nextPointIndex], isSoftBody?stiffness:500));
            springs.push_back(Spring(points[pointIndex], points[nextNextPointIndex], isSoftBody?stiffness:500));
        }
    }

private:
    void CreateEdges()
    {
        for (size_t i = 0; i < points.size(); i++)
        {
            int pointIndex = i;
            int nextPointIndex = (i + 1) % points.size();

            edges.push_back(Edge(&points[pointIndex], &points[nextPointIndex]));
        }
    }

public:
    void ApplyTorque(float t)
    {
        torque += t;
    }

public:
    void ApplyForce(Vector2f t)
    {
        addForce += t;
    }
};


// Function to create a thick line as a rectangle
RectangleShape CreateThickLine(Vector2f start, Vector2f end, float thickness, Color color) {
    Vector2f direction = end - start;
    float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
    
    RectangleShape line(Vector2f(length, thickness));
    line.setFillColor(color);
    line.setPosition(start);
    
    // Calculate angle and rotate
    float angleRadians = std::atan2(direction.y, direction.x);
    line.setRotation(sf::radians(angleRadians));
    
    //// Offset to center the thickness
    //sf::Vector2f offset(-thickness/2.0f * std::sin(angleRadians), 
    //                    thickness/2.0f * std::cos(angleRadians));
    //line.setPosition(start + offset);
    
    return line;
}









void DrawShapeOnScreen(RenderWindow& window, const Shape& shape) 
{
    for (size_t i = 0; i < shape.points.size(); i++)
    {
        int pointIndex = i;
        int nextPointIndex = (i + 1) % shape.points.size();
        int nextNextPointIndex = (i + 2) % shape.points.size();


        RectangleShape line = CreateThickLine(shape.points[pointIndex].pos, shape.points[nextPointIndex].pos, 5.0, shape.edges[i].color);
        RectangleShape line2 = CreateThickLine(shape.points[pointIndex].pos, shape.points[nextNextPointIndex].pos, .1f, Color::Yellow);

        window.draw(line);
        window.draw(line2);
    }


    //---------Create Vertices----------//

    CircleShape pointCircle;
    pointCircle.setRadius(8.0f);
    pointCircle.setFillColor(Color::White);
    pointCircle.setOrigin({ 8.0f, 8.0f });

    for (const auto& point : shape.points) {
        pointCircle.setPosition(point.pos);
        pointCircle.setFillColor(point.color);
        window.draw(pointCircle);
    }
}



bool IsPointInRectangle(const Vector2f& point, const Shape& rect) {
    // Assuming rectangle points are in order: top-left, top-right, bottom-right, bottom-left
    float minX = std::min({rect.points[0].pos.x, rect.points[1].pos.x, rect.points[2].pos.x, rect.points[3].pos.x});
    float maxX = std::max({rect.points[0].pos.x, rect.points[1].pos.x, rect.points[2].pos.x, rect.points[3].pos.x});
    float minY = std::min({rect.points[0].pos.y, rect.points[1].pos.y, rect.points[2].pos.y, rect.points[3].pos.y});
    float maxY = std::max({rect.points[0].pos.y, rect.points[1].pos.y, rect.points[2].pos.y, rect.points[3].pos.y});
    
    return (point.x >= minX && point.x <= maxX && point.y >= minY && point.y <= maxY);
}

void HandleShapeCollision(Shape& innerShape, const Shape& outerShape) {
    float minX = std::min({ outerShape.points[0].pos.x, outerShape.points[1].pos.x, outerShape.points[2].pos.x, outerShape.points[3].pos.x });
    float maxX = std::max({ outerShape.points[0].pos.x, outerShape.points[1].pos.x, outerShape.points[2].pos.x, outerShape.points[3].pos.x });
    float minY = std::min({ outerShape.points[0].pos.y, outerShape.points[1].pos.y, outerShape.points[2].pos.y, outerShape.points[3].pos.y });
    float maxY = std::max({ outerShape.points[0].pos.y, outerShape.points[1].pos.y, outerShape.points[2].pos.y, outerShape.points[3].pos.y });

    for (auto& point : innerShape.points) {
        bool collided = false;

        if (point.pos.x < minX) { point.pos.x = minX; point.velocity.x *= -0.5f; collided = true; }
        if (point.pos.x > maxX) { point.pos.x = maxX; point.velocity.x *= -0.5f; collided = true; }
        if (point.pos.y < minY) { point.pos.y = minY; point.velocity.y *= -0.5f; collided = true; }
        if (point.pos.y > maxY) { point.pos.y = maxY; point.velocity.y *= -0.5f; collided = true; }

        if (collided) {
            // add torque effect from bounce impulse
            Vector2f r = point.pos - innerShape.centerOfMass;
            Vector2f impulse = point.velocity * point.mass;
            float torque = r.x * impulse.y - r.y * impulse.x;
            innerShape.torque += torque;
        }
    }
}

struct ClosestEdgeToPoint
{
    Edge* edge = nullptr;
    float distance = FLT_MAX;
    bool inside = false;
    Vector2f directionToEdge = Vector2f();
	Vector2f closestPoint = Vector2f();
    Vector2f edgeDir = Vector2f();
};

ClosestEdgeToPoint GetClosestEdge(const Vector2f& point, Shape& shape) 
{
    ClosestEdgeToPoint result;
    bool hasEdge = false;
    
    // First pass: find the closest edge and determine if point is inside
    int intersections = 0;
    
    for (Edge& edge : shape.edges)
    {
        Vector2f p1 = edge.point1->pos;
        Vector2f p2 = edge.point2->pos;
        
        // Check for ray intersection (for point-in-polygon test)
        if ((p1.y > point.y) != (p2.y > point.y)) {
            float x_intersect = (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
            if (x_intersect > point.x) {
                intersections++;
            }
        }
        
        // Find closest point on this edge
        Vector2f edgeVec = p2 - p1;
        float edgeLength = edgeVec.length();
        
        if (edgeLength < 0.0001f) continue; // Skip zero-length edges
        
        Vector2f edgeDir = edgeVec / edgeLength;
        Vector2f p1toPoint = point - p1;
        
        // Project point onto edge
        float projection = p1toPoint.dot(edgeDir);
        projection = std::max(0.0f, std::min(edgeLength, projection));
        
        // Find closest point on edge
        Vector2f closestPoint = p1 + edgeDir * projection;
        Vector2f toPoint = closestPoint - point;
        float distance = toPoint.length();
        
        if (!hasEdge || distance < result.distance) {
            // Calculate direction to edge (pointing away from the shape)
            Vector2f direction = closestPoint - point;
            float dirLength = direction.length();
            
            if (dirLength > 0.0001f) {
                result.directionToEdge = direction / dirLength;
            } else {
                // If point is on the edge, use edge normal
                result.directionToEdge = Vector2f(-edgeDir.y, edgeDir.x);
            }
            
            result.distance = distance;
            result.edge = &edge;
            result.closestPoint = closestPoint;
            result.edgeDir = edgeDir;
            hasEdge = true;
        }else if (!hasEdge || distance == result.distance)
        {
            Vector2f direction = closestPoint - point;
            float dirLength = direction.length();
            if (abs(direction.dot(edgeDir)) < abs(direction.dot(result.edgeDir)))
            {
                if (dirLength > 0.0001f) {
                    result.directionToEdge = direction / dirLength;
                } else {
                    // If point is on the edge, use edge normal
                    result.directionToEdge = Vector2f(-edgeDir.y, edgeDir.x);
                }
                
                result.distance = distance;
                result.edge = &edge;
                result.closestPoint = closestPoint;
                result.edgeDir = edgeDir;
                hasEdge = true;
            }
        }
    }
    
    // Point is inside if number of intersections is odd
    result.inside = (intersections % 2 == 1);
    
    return result;
}

void CheckForCollision(Vector2f point, Edge* edge)
{
    Vector2f p1 = edge->point1->pos;
    Vector2f p2 = edge->point2->pos;
    
    Vector2f edgeVec = p2 - p1;
        float edgeLength = edgeVec.length();
        
        if (edgeLength < 0.0001f) return; // Skip zero-length edges
        
        Vector2f edgeDir = edgeVec / edgeLength;
        Vector2f p1toPoint = point - p1;
        
        // Project point onto edge
        float projection = p1toPoint.dot(edgeDir);
        projection = std::max(0.0f, std::min(edgeLength, projection));
        
        // Find closest point on edge
        Vector2f closestPoint = p1 + edgeDir * projection;
        Vector2f toPoint = closestPoint - point;
        float distance = toPoint.length();
        
        if (distance < 15.0f)
        {
            std::cout << "Collision detected" << std::endl;   
        }
}


void HandleSoftBodyCollision(Shape& shapeA, Shape& shapeB)
{
    const float collisionDistance = 5.0f;
    const float restitution = 1.f; 
    const float friction = 0.4f;    
    int i = 0;
    for (auto& pointA : shapeA.points) {
        // Skip if point has infinite mass (static)
        if (pointA.mass <= 0.0f) continue;
        
        // Get closest edge on shapeB to this point
        ClosestEdgeToPoint closestEdge = GetClosestEdge(pointA.pos, shapeB);



        if (!closestEdge.edge) continue; // Skip if no valid edge found

        // Calculate edge properties
        Vector2f edgeVec = closestEdge.edge->point2->pos - closestEdge.edge->point1->pos;
        float edgeLength = edgeVec.length();
        if (edgeLength < 0.0001f) continue; // Skip zero-length edges
        
        Vector2f edgeDir = closestEdge.edgeDir;
        Vector2f directionToEdge = closestEdge.directionToEdge;
        
        // Calculate projection of point onto edge
        Vector2f p1ToPoint = pointA.pos - closestEdge.edge->point1->pos;
        float projection = p1ToPoint.dot(edgeDir);
        
        // Calculate collision depth
        float depth = collisionDistance - abs(closestEdge.distance);


        if (closestEdge.inside)
        {
            pointA.color = Color::Red;


            // Take it outt
            Vector2f correctionVec = directionToEdge * (abs(closestEdge.distance) + collisionDistance * 0.5f)     *     0.2f;
            Vector2f correctionVec = directionToEdge * (abs(closestEdge.distance) + collisionDistance * 0.5f)    *    0.2f;
            
            pointA.pos += correctionVec / 3.f;
			closestEdge.edge->point1->pos -= correctionVec / 3.f * 2.f;
			closestEdge.edge->point2->pos -= correctionVec / 3.f * 2.f;
            //CheckForCollision(pointA.pos, closestEdge.edge);

            continue;
        }else
        {
            pointA.color = Color::White;
        }

        //visual
        if (i == 0) {
            for (auto& edge : shapeB.edges)
            {
                edge.color = Color::White;
            }
            closestEdge.edge->color = Color::Cyan;
            pointA.color = Color::Cyan;
        }
        i++;

        // Skip if projection is outside edge bounds
        if (projection < 0 || projection > edgeLength) continue;
        
        if (collisionDistance < closestEdge.distance)
        {
            continue;
        }

        // Skip if edge direction is invalid
        if (std::isnan(projection)) {
            continue;
        }
        
        float e = std::min(restitution, 1.0f);
        
        Vector2f pointVelocity = pointA.velocity;
        Vector2f pointVelocityPerpen = pointVelocity.dot(directionToEdge) * directionToEdge;
        float pointMass = pointA.mass;
        //
        float mass1 = closestEdge.edge->point1->mass;
        float mass2 = closestEdge.edge->point2->mass;
        Vector2f velo1 = closestEdge.edge->point1->velocity;
        Vector2f velo2 = closestEdge.edge->point2->velocity;
        Vector2f edgeVelocity = (velo1 + velo2) / 2.f;
        Vector2f edgeVelocityPerpen = edgeVelocity.dot(directionToEdge) * directionToEdge;

        Vector2f relativeVel = pointVelocity - edgeVelocity;
        Vector2f relativeVelPerpen = relativeVel.dot(directionToEdge) * directionToEdge;

        Vector2f subVelPoint = (2.f * relativeVelPerpen / 3.f) * (1 + e);
        Vector2f addVelEdge = (relativeVelPerpen / 3.f) * (1 + e);

        if (pointVelocity.dot(directionToEdge) < -0.1f) continue; // Small threshold to prevent jitter


        float ratio = projection / edgeLength;

        pointA.velocity -= subVelPoint;
        closestEdge.edge->point1->velocity += 2.f * addVelEdge * (1 - ratio);
        closestEdge.edge->point2->velocity += 2.f * addVelEdge * ratio;

        //Apply torque
        Vector2f r = pointA.pos - shapeA.centerOfMass;
        float torque = - r.x * subVelPoint.y + r.y * subVelPoint.x;
        shapeA.torque += torque; 

        //Apply torque to edge
        r = closestEdge.edge->point1->pos - shapeA.centerOfMass;
        torque = r.x * addVelEdge.y - r.y * addVelEdge.x;
        shapeA.torque += torque; 

        r = closestEdge.edge->point2->pos - shapeA.centerOfMass;
        torque = r.x * addVelEdge.y - r.y * addVelEdge.x;
        shapeA.torque += torque; 
    }
}



void CalculatePhysics(Shape& shape)
{
    for (auto& point : shape.points)
        point.force = { 0.f, 0.f };

    shape.CalculateCenterOfMass();


    for (auto& spring : shape.springs)
    {
        Vector2f delta = spring.point2.pos - spring.point1.pos;
        float dist = std::sqrt(delta.x * delta.x + delta.y * delta.y);

        if (dist > 0.0001f)
        {
            Vector2f dir = delta / dist; // normalize
            float displacement = dist - spring.restLength;

            // Hooke�s law
            Vector2f force = dir * (spring.stiffnes * displacement);

            // Apply equal and opposite forces
            spring.point1.force += force;
            spring.point2.force -= force;
        }
    }

    for (auto& point : shape.points)
    {
        if (shape.staticBody) continue;

        point.force += {0, 9.8f};
        point.force += point.addForce;
        point.force += shape.addForce / (float)shape.points.size();
        point.addForce = {0, 0};

        if (point.force != Vector2f(0,0))
        {
            point.acceleration = point.force / point.mass;
        }

        point.velocity += (point.acceleration * 1.0f / 100.0f );
        if (point.velocity.length() < 0.0001f) point.velocity = { 0,0 };


        point.pos += point.velocity * 1.0f / 100.0f;

        point.velocity *= 0.999f;   


        //-----Add Torque From Force-----//
        Vector2f r = point.pos - shape.centerOfMass;
        float torque = r.x * point.force.y - r.y * point.force.x; 
        shape.torque += torque; 
    }
    shape.addForce = { 0,0 };

    //--------NOw Torque------//
    if (!shape.staticBody && shape.momentOfInertia > 0.f)
    {
        float angularAcc = shape.torque / shape.momentOfInertia;
        shape.angularVelocity += angularAcc * 1.0f / 100.0f;
        shape.angularVelocity *= 0.999f;

        float angle = shape.angularVelocity * 1.0f / 100.0f;

        //------------Rotate all points around center------//
        if (std::abs(angle) > 0.0001) {
            for (auto& point : shape.points)
            {
                Vector2f r = point.pos - shape.centerOfMass;
                float initialAngle = atan2f(r.y, r.x);
                float finalAngle = initialAngle + angle;

                Vector2f rotated(r.length() * cos(finalAngle), r.length() * sin(finalAngle));

                point.pos = shape.centerOfMass + rotated;
            }
        }
        shape.torque = 0.f; // reset torque each frame
    }
}




int main() {
    //======================Create Window=======================//
    RenderWindow window(VideoMode({ 1920, 1080 }), "Soft Smooth Body ;)");





    //======================Initialize========================//
    
    
    //======================FPS Tracker=========================//
    using clock = std::chrono::high_resolution_clock;
    auto lastTime = clock::now();
    int frames = 0;
    

    vector<Point> points;

    points.push_back({ {500.f, 300.f} }); // top
    //points.push_back({ {579.f, 317.f} });
    //points.push_back({ {641.f, 370.f} });
    points.push_back({ {700.f, 300.f} });
    //points.push_back({ {679.f, 450.f} });
    points.push_back({ {700.f, 500.f} });
    //points.push_back({ {641.f, 630.f} });
    //points.push_back({ {579.f, 683.f} });
    points.push_back({ {500.f, 500.f} }); // bottom
    //points.push_back({ {421.f, 683.f} });
    //points.push_back({ {359.f, 630.f} });
    points.push_back({ {400.f, 400.f} });
    //points.push_back({ {359.f, 370.f} });
    //points.push_back({ {421.f, 317.f} });
    //points.push_back({ {500.f, 300.f} }); // repeat first point if you want closed loop
    Shape shape(points, true);

    // ----- second shape placed next to first (shifted by +250 on X axis) -----
    Point q1 = { {350.f, 100.f} };
    Point q2 = { {350.f, 200.f} };
    Point q3 = { {800.f, 180.f} };
    Point q4 = { {750.f, 200.f} };
    Point q5 = { {750.f, 100.f} };
    std::vector<Point> points2 = { q1, q2, /*q3, */q4, q5 };
    Shape shapeB(points2, true);  // second shape
    Shape shapeC(points2, true);  // second shape

    // Beveled square centered at (400,150), width=100, bevel=20
// Adjust coordinates to match your setup
    Point b1 = { {350.f, 100.f} }; // bottom-left bevel start
    Point b2 = { {370.f, 100.f} };
    Point b3 = { {450.f, 100.f} }; // bottom-right bevel start
    Point b4 = { {450.f, 120.f} };
    Point b5 = { {450.f, 200.f} }; // top-right bevel start
    Point b6 = { {430.f, 200.f} };
    Point b7 = { {350.f, 200.f} }; // top-left bevel start
    Point b8 = { {350.f, 180.f} };

    points2 = { b1, b2, b3, b4, b5, b6, b7, b8 };
    //Shape shapeC(points2, false);  // beveled square shape



    Point p1 = { {50.f, 50.f} };
    Point p2 = { {1500.f, 50.f} };
    Point p3 = { {1500.f, 1000.f} };
    Point p4 = { {50.f, 1000.f} };

    points.clear();
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);

    Shape boundaryShape(points, false);

    boundaryShape.staticBody = true;






    
    //======================Game Loop=======================//

    while (window.isOpen()) {
        // Event handling
        while (auto event = window.pollEvent()) {
            // event is an std::optional<sf::Event>
            const Event ev = *event;

            // Check if user closed the window
            if (ev.is<Event::Closed>()) {
                window.close();
            }
        }

        CalculatePhysics(shape);
        CalculatePhysics(shapeB);
        CalculatePhysics(shapeC);

        HandleShapeCollision(shape, boundaryShape);
        HandleShapeCollision(shapeB, boundaryShape);
        HandleShapeCollision(shapeC, boundaryShape);

		HandleSoftBodyCollision(shape, shapeB);
		HandleSoftBodyCollision(shapeB, shape);
		HandleSoftBodyCollision(shapeC, shape);
		HandleSoftBodyCollision(shape, shapeC);
		HandleSoftBodyCollision(shapeC, shapeB);
		HandleSoftBodyCollision(shapeB, shapeC);



        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Up)) {
            shape.ApplyTorque(-10000.f);
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Down)) {
            shape.ApplyTorque(10000.f);
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) {
            shape.ApplyForce({ -50, 0 });
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) {
            shape.ApplyForce({ 50, 0 });
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space)) {
            shape.ApplyForce({0, -300.f });
        }






        

        //===================Draw Visual=======================//
        window.clear(Color::Black);
        //DrawShapeWithPoints(window, shape);
        DrawShapeOnScreen(window, shape);
        DrawShapeOnScreen(window, boundaryShape);
        DrawShapeOnScreen(window, shapeB);
        DrawShapeOnScreen(window, shapeC);

        
        window.setFramerateLimit(1000);
        window.display();



        // FPS counter
        frames++;
        auto now = clock::now();
        std::chrono::duration<float> elapsed = now - lastTime;

        if (elapsed.count() >= 1.0f) { // every second
            std::cout << "FPS: " << frames << std::endl;
            frames = 0;
            lastTime = now;
        }
    }









    //======================End Game Loop=======================//

    return 0;
}




// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
