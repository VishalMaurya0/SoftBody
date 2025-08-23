// FirstC++.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <SFML/Graphics.hpp>
#include <optional> // needed for std::optional
#include <cmath> // needed for sqrt, atan2, sin, cos
#include "SoftBody.h"

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



struct Shape {
    vector<Point> points;
    Vector2f centerOfMass;
    Color color = Color::White;

    bool isSoftBody = true;
    bool staticBody = false;

    vector<Spring> springs;

    float angularVelocity = 0.f;
    float torque = 0.f;
    float momentOfInertia = 0.f;


    Shape(const vector<Point>& points, const bool isSoftBody) {
        this->isSoftBody = isSoftBody;
        this->points = points;
        CalculateCenterOfMass();
        CreateSprings();
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

private:
    void CreateSprings()
    {
        for (size_t i = 0; i < points.size(); i++)
        {
            int pointIndex = i;
            int nextPointIndex = (i + 1) % points.size();
            int nextNextPointIndex = (i + 2) % points.size();

            springs.push_back(Spring(points[pointIndex], points[nextPointIndex], isSoftBody?3:500));
            springs.push_back(Spring(points[pointIndex], points[nextNextPointIndex], isSoftBody?3:500));
        }
    }

public:
    void ApplyTorque(float t)
    {
        torque += t;
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


        RectangleShape line = CreateThickLine(shape.points[pointIndex].pos, shape.points[nextPointIndex].pos, 5.0, Color::White);
        RectangleShape line2 = CreateThickLine(shape.points[pointIndex].pos, shape.points[nextNextPointIndex].pos, 2.0, Color::Yellow);

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
    Point p1 = { {500.f, 350.f} };   // top
    Point p2 = { {400.f, 400.f} };   // top-left
    Point p3 = { {400.f, 600.f} };   // bottom-left
    Point p4 = { {500.f, 650.f} };   // bottom
    Point p5 = { {600.f, 600.f} };   // bottom-right
    Point p6 = { {600.f, 400.f} };   // top-right

    vector<Point> points;

    //points.push_back({ {500.f, 300.f} }); // top
    //points.push_back({ {579.f, 317.f} });
    //points.push_back({ {641.f, 370.f} });
    points.push_back({ {679.f, 450.f} });
    points.push_back({ {679.f, 550.f} });
    points.push_back({ {641.f, 630.f} });
    points.push_back({ {579.f, 683.f} });
    //points.push_back({ {500.f, 700.f} }); // bottom
    //points.push_back({ {421.f, 683.f} });
    //points.push_back({ {359.f, 630.f} });
    points.push_back({ {321.f, 550.f} });
    points.push_back({ {321.f, 450.f} });
    points.push_back({ {359.f, 370.f} });
    points.push_back({ {421.f, 317.f} });
    //points.push_back({ {500.f, 300.f} }); // repeat first point if you want closed loop


    Shape shape(points, true);

    p1 = { {50.f, 50.f} };
    p2 = { {1000.f, 50.f} };
    p3 = { {1000.f, 1000.f} };
    p4 = { {50.f, 1000.f} };

    points.clear();
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);

    Shape shape2(points, false);

    shape2.staticBody = true;






    
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
        HandleShapeCollision(shape, shape2);

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left)) {
            shape.ApplyTorque(-2000.f);
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Right)) {
            shape.ApplyTorque(2000.f);
        }






        

        //===================Draw Visual=======================//
        window.clear(Color::Black);
        //DrawShapeWithPoints(window, shape);
        DrawShapeOnScreen(window, shape);
        DrawShapeOnScreen(window, shape2);
        window.display();
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
