#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <cfloat>
#include <tuple>
#include <omp.h>
#include <list>
#include <chrono>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


/*  ------------------------ VECTOR ---------------------------------   */

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    Vector normalize() {
        double n = norm();
        return Vector(data[0]/n, data[1]/n, data[2]/n);
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
    public:
        Vector O;
        Vector u;
        Ray(const Vector& O, const Vector& u) : O(O), u(u) {};

};

struct Intersection {
    Vector P;
    Vector N;
    Vector albedo = Vector(0., 0., 0.);
    double t = 0;
    double refraction_index = 1;
    bool reflective = false;
    bool intersects = false;
    Intersection(Vector albedo = Vector(0., 0., 0.), bool reflective = false, double refraction_index = 1) : albedo(albedo), refraction_index(refraction_index), reflective(reflective) {};
};

class Geometry {
public:
    virtual Intersection intersect(const Ray& ray) = 0;
};

class Sphere : public Geometry {
    private:
        Vector C;
        Vector albedo;
        double R;
        double refraction_index;
        bool reflective;
        bool hollow;

    public:
        explicit Sphere(Vector C, double R, Vector albedo, bool reflective = false, double refraction_index = 1., bool hollow = false) : C(C), R(R), albedo(albedo), reflective(reflective), refraction_index(refraction_index), hollow(hollow) {};

        Intersection intersect(const Ray& r) override {
            Intersection i(this->albedo, this->reflective, this->refraction_index);

            Vector origin_center = r.O - this->C;
            double delta=dot(r.u, r.O-C)*dot(r.u, r.O-C) - ((r.O-C).norm2()-R*R);
            if (delta >= 0.) {
                double t1=dot(r.u, C-r.O)-sqrt(delta);
				double t2=dot(r.u, C-r.O)+sqrt(delta);
                if(t2>=0 && t1>=0){
					i.t=t1;
                    i.intersects=true;
				}
				else if (t2>=0 && t1<0){
					i.t=t2;
                    i.intersects=true;
				}
                else{
                    i.t=0;
                    i.intersects=false;
                }
            }

            i.P = r.O + r.u * i.t;
            i.N = (i.P - C).normalize();
            i.N = (!hollow)? i.N : (-1)*i.N;
            return i;
        }
};

class BoundingBox {
public:
    Vector B_min;
    Vector B_max;

    explicit BoundingBox(Vector m = Vector(), Vector M = Vector()) {
        B_min = m;
        B_max = M;
    }
};

struct Node {
    BoundingBox bounding_box;
    int first_triangle;
    int last_triangle;
    Node* left_child;
    Node* right_child;
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};


class TriangleMesh : public Geometry {
    double sf;
    Vector translation;
    Vector albedo;
    double refraction_index;
    bool reflective;

public:

    ~TriangleMesh() {}
    TriangleMesh(double sf, Vector translation, Vector albedo = Vector(0, 0, 0), double refraction_index = 1, bool reflective = false) : sf(sf), translation(translation), albedo(albedo), refraction_index(refraction_index), reflective(reflective), root(new Node) {};

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);
        this->bvh_tree(root, 0, indices.size());
    }

    BoundingBox compute_bounding_box(int first_triangle, int last_triangle) {
        double mX = 100000000, MX = -100000000,  MY = -100000000, mY = 100000000, mZ = 100000000, MZ = -100000000;

        for (int i=first_triangle; i<last_triangle; i++) {
            auto triangle_vertices = {vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk]};
            for (auto& vertex : triangle_vertices) {
                Vector scaled_translated_vertex = sf * vertex + translation;
                mX = std::min(mX, scaled_translated_vertex[0]);
                mY = std::min(mY, scaled_translated_vertex[1]);
                mZ = std::min(mZ, scaled_translated_vertex[2]);

                MX = std::max(MX, scaled_translated_vertex[0]);
                MY = std::max(MY, scaled_translated_vertex[1]);
                MZ = std::max(MZ, scaled_translated_vertex[2]);
            }
        }
        return BoundingBox(Vector(mX, mY, mZ), Vector(MX, MY, MZ));
    }

    bool intersect_bounding_box(const Ray& ray, BoundingBox bounding_box, double* t) {
        double tx1, ty1, tz1, tx0, ty0, tz0, t_B_min, t_B_max;
        Vector V1 = Vector(1, 0, 0), V2 = Vector(0, 1, 0), V3 = Vector(0, 0, 1);
        
        t_B_min = dot(bounding_box.B_min - ray.O, V1) / dot(ray.u, V1);
        t_B_max = dot(bounding_box.B_max - ray.O, V1) / dot(ray.u, V1);
        tx0 = std::min(t_B_min, t_B_max);
        tx1 = std::max(t_B_min, t_B_max);

        t_B_min = dot(bounding_box.B_min - ray.O, V2) / dot(ray.u, V2);
        t_B_max = dot(bounding_box.B_max - ray.O, V2) / dot(ray.u, V2);
        ty0 = std::min(t_B_min, t_B_max);
        ty1 = std::max(t_B_min, t_B_max);

        t_B_min = dot(bounding_box.B_min - ray.O, V3) / dot(ray.u, V3);
        t_B_max = dot(bounding_box.B_max - ray.O, V3) / dot(ray.u, V3);
        tz0 = std::min(t_B_min, t_B_max);
        tz1 = std::max(t_B_min, t_B_max);

        auto m = std::min(std::min(tx1,ty1),tz1), M = std::max(std::max(tx0,ty0),tz0);
        if(m > M > 0){
            *t = M;
            return true;
        }
        return false;
    }

    void bvh_tree(Node *node, int first_triangle, int last_triangle) {
        node->bounding_box = this->compute_bounding_box(first_triangle, last_triangle);
        node->first_triangle = first_triangle;
        node->last_triangle = last_triangle;

        Vector diag = node->bounding_box.B_max - node->bounding_box.B_min;
        int longest_axis=0;
        auto M = abs(diag[0]);
        for(int i = 1; i < 3; i++){
            if(abs(diag[i]) > M){
                M = abs(diag[i]);
                longest_axis = i;
            }
        }

        Vector middle = node->bounding_box.B_min + diag*0.5;
        int pivot = first_triangle;
        Vector v1, v2, v3, BaryCenter;
        for (int i = first_triangle; i < last_triangle; i++) {
            v1 = (sf * vertices[indices[i].vtxi]) + translation;
            v2 = (sf * vertices[indices[i].vtxj]) + translation;
            v3 = (sf * vertices[indices[i].vtxk]) + translation;
            BaryCenter = (v1 + v2 + v3)/3;
            if(BaryCenter[longest_axis] < middle[longest_axis]){
                std::swap(indices[i], indices[pivot]);
                pivot++;
            }
        }

        if (pivot <= first_triangle || pivot >= last_triangle - 1){
            return;
        }

        node->left_child = new Node();
        node->right_child = new Node();
        bvh_tree(node->left_child, first_triangle, pivot);
        bvh_tree(node->right_child, pivot, last_triangle);
    }

    Intersection intersect(const Ray &ray) override{

        Intersection i(this->albedo, this->reflective, this->refraction_index);
        double t_min = 1000000000, t;
        if(!intersect_bounding_box(ray, this->root->bounding_box, &t)){
            return Intersection();
        }

        //no BVH -------------------------------
        // Vector A, B, C, N; 
        // double Nu, alpha, beta, gamma;
        // for (int j = 0; j < indices.size(); j++){
        //     TriangleIndices triangle_indices = indices[j];
        //     A = sf * vertices[triangle_indices.vtxi] + translation;
        //     B = sf * vertices[triangle_indices.vtxj] + translation;
        //     C = sf * vertices[triangle_indices.vtxk] + translation;
        //     N = cross(B - A, C - A);
        
        //     Nu = dot(ray.u, N);
        //     beta = dot(C - A, cross(A - ray.O, ray.u)) / Nu;
        //     gamma = - dot(B - A, cross(A - ray.O, ray.u)) / Nu;
        //     alpha = 1 - beta - gamma;
        //     t = dot(A - ray.O, N) / Nu;

        //     if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < t_min) {
        //         t_min = t;
        //         i.intersects = true;
        //         i.t = t;
        //         i.P = A + beta * (B - A) + gamma * (C - A);
        //         i.N = N;
        //     }
        // }
        // return i;
        /*----------------------------------*/
        //With BVH
        std::list<Node*> nodes;
        nodes.push_front(root);
        while(!nodes.empty()){
            Node *currNode = nodes.back();
            nodes.pop_back();

            if(currNode->left_child){
                if(intersect_bounding_box(ray, currNode->left_child->bounding_box, &t) && t < t_min){
                        nodes.push_back(currNode->left_child);
                }
                if(intersect_bounding_box(ray, currNode->right_child->bounding_box, &t) && t < t_min){
                        nodes.push_back(currNode->right_child);
                }
            } 
            else{
                Vector A, B, C, N; 
                double Nu, alpha, beta, gamma;
                for (int j = currNode->first_triangle; j < currNode->last_triangle; j++){
                    TriangleIndices triangle_indices = indices[j];
                    A = sf * vertices[triangle_indices.vtxi] + translation;
                    B = sf * vertices[triangle_indices.vtxj] + translation;
                    C = sf * vertices[triangle_indices.vtxk] + translation;
                    N = cross(B - A, C - A);
                
                    Nu = dot(ray.u, N);
                    beta = dot(C - A, cross(A - ray.O, ray.u)) / Nu;
                    gamma = - dot(B - A, cross(A - ray.O, ray.u)) / Nu;
                    alpha = 1 - beta - gamma;
                    t = dot(A - ray.O, N) / Nu;

                    if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < t_min) {
                        t_min = t;
                        i.intersects = true;
                        i.t = t;
                        i.P = A + beta * (B - A) + gamma * (C - A);
                        i.N = N;
                    }
                }
            }
        }
        return i;
    }
    
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox bounding_box;
    Node* root;
};


class Scene {
    private:
        std::vector<Geometry*> geometries;
        Vector source;
        double intensity = 1e5;

    public:
        explicit Scene(Vector source) : source(source){};

        void addGeometry(Geometry* geometry){ 
            geometries.push_back(geometry); 
        }

        //Auxiliary random_cos function as per the formulas in lecture 2
        Vector random_cos(const Vector &N){
            double m = 10000000;
            int index = 0;
            for(int i = 0; i < 3; i++){
                if(abs(N[i]) < m){
                    m = abs(N[i]);
                    index = i;
                }
            }

            Vector T1, T2;
            if(index==0){
                T1 = Vector(0, N[2], -N[1]).normalize();
            }
            else if (index==1){
                T1 = Vector(N[2], 0., -N[0]).normalize();
            }
            else{
                T1 = Vector(N[1], -N[0], 0.).normalize();
            }
            T2 = cross(N, T1);

            double r1 = ((double) rand() / (RAND_MAX));
            double r2 = ((double) rand() / (RAND_MAX));

            return T1 * sqrt(1-r2) * cos(2 * M_PI * r1) + T2 * sqrt(1-r2) * sin(2 * M_PI * r1) + N * sqrt(r2);
        }

        //Find the correct intersection
        Intersection intersect(const Ray& ray){
            Intersection first_intersection, found_intersection;
            double t = 1000000000;
            for(auto& geometry : geometries){
                found_intersection = geometry->intersect(ray);
                if(found_intersection.intersects && found_intersection.t < t){
                    t = found_intersection.t;
                    first_intersection = found_intersection;
                }
            }
            return first_intersection;
        }


        Vector getColor(const Ray& ray, int depth) {
            if(depth < 0){
                return Vector(0., 0., 0.);
            }
            Intersection i = intersect(ray);
            Vector Lo(0, 0, 0);
            if(i.intersects){
                //handle reflection (lec. 1 slide 37)
                Vector P = i.P + (i.N * 1e-10);
                if (i.reflective) {
                    Ray reflected_ray = Ray(P, ray.u - (2 * dot(ray.u, i.N) * i.N));
                    return getColor(reflected_ray, depth - 1);
                }

                //handle refraction (lec.1 slides 39-14)
                if(i.refraction_index != 1){
                    double Nu = dot(ray.u, i.N);

                    double n1,n2;
                    if(Nu >= 0){
                        n1=i.refraction_index;
                    }
                    else{
                        n1=1;
                    }
                    if(Nu >= 0){
                        n2=1;
                    }
                    else{
                        n2 = i.refraction_index;
                    }
                    Vector invertN = (Nu > 0) ? -1.*i.N : i.N;

                    P = i.P - (invertN * 1e-10);
                    Nu = dot(ray.u, invertN);
                    
                    //Fresnel
                    if(1 - (n1 * n1) / (n2 * n2) * (1 - Nu * Nu) > 0){
                        Vector u_n = -1. * invertN * sqrt(1 - (n1 * n1) / (n2 * n2) * (1 - Nu * Nu));
                        Vector u_t = (n1/n2) * (ray.u - Nu * invertN);
                        Vector u = u_t + u_n;
                        double R = pow((n1-n2)/(n1+n2), 2) + (1-pow((n1-n2)/(n1+n2), 2)) * pow(1 - abs(dot(invertN, u)), 5);
                        if ((double) rand() / (RAND_MAX) < R) {
                            Ray reflected_ray = Ray(P, ray.u - (2 * i.N * Nu));
                            return getColor(reflected_ray, depth - 1);
                        } else {
                            Ray refracted_ray = Ray(P, u);
                            return getColor(refracted_ray, depth - 1);
                        }
                    } 
                    else{
                        Ray reflected_ray = Ray(P, ray.u - (2 * dot(i.N, ray.u) * i.N));
                        return getColor(reflected_ray, depth - 1);
                    }
                }
                double distance = (source - P).norm();
                Vector u_i = (source - P).normalize();
                Intersection i_light = intersect(Ray(source, u_i*(-1)));
                bool visible; 
                if(!i_light.intersects || i_light.t > distance){
                    visible = 1;
                } 
                else{
                    visible=0;
                }
                Lo = intensity / (4*M_PI*distance*distance) * i.albedo / M_PI * (double)visible * std::max((double)0, dot(u_i, i.N));

                //indirect lighting
                Ray random_ray = Ray(P, random_cos(i.N));
                Lo = Lo + i.albedo*getColor(random_ray, depth - 1);
            }
            return Lo;
        }
};



void BoxMuller(double sigma, double& x, double& y) {
    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));
    x = sigma * sqrt(-2*log(r1)) * cos(2*M_PI*r2);
    y = sigma * sqrt(-2*log(r1)) * sin(2*M_PI*r2);
}

int main() {

    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    Scene scene = Scene(Vector(-10, 20, 40));

    // reflective, refractive and hollow spheres
    // Sphere reflective = Sphere(Vector(-20, 0, 0), 10, Vector(1., 1., 1.), true);
    // Sphere refractive = Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    // Sphere exterior_hollow = Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    // Sphere interior_hollow = Sphere(Vector(20, 0, 0), 9.99, Vector(1., 1., 1.), false, 1.5, true);


    //floor, ceiling and walls created as massive spheres, as suggested in the lecture
    Sphere ceiling = Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere floor = Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere front = Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere back = Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere left = Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
    Sphere right = Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));

    //add walls and spheres to the scene
    // scene.addGeometry(&reflective);
    // scene.addGeometry(&refractive);
    // scene.addGeometry(&exterior_hollow);
    // scene.addGeometry(&interior_hollow);
    scene.addGeometry(&ceiling);
    scene.addGeometry(&floor);
    scene.addGeometry(&front);
    scene.addGeometry(&back);
    scene.addGeometry(&left);
    scene.addGeometry(&right);

    // add cat to the scene
    TriangleMesh cat = TriangleMesh(0.3, Vector(0, -10, 20), Vector(1., 1., 1.));
    cat.readOBJ("cat.obj");
    scene.addGeometry(&cat);

    int W = 512;
    int H = 512;
    std::vector<unsigned char> image(W*H*3, 0);
    Vector Camera = Vector(0, 0, 55);
    double fov = 1.0472;
    double gamma = 2.2;
    int max_depth = 5;
    int rays_per_pixel = 10;

    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < H; i++){
        for(int j = 0; j < W; j++){
            Vector pixelColor = Vector(0, 0, 0);
            double x, y;
            for (int k = 0; k < rays_per_pixel; k++){
                BoxMuller(0.5, x, y);
                //formulas on page 15 lecture notes
                Vector Pixel = Vector(Camera[0]+(j+x)+0.5-W/2, Camera[1]-(i+y)-0.5+H/2, Camera[2]-W/(2*tan(fov/2)));
                Ray ray = Ray(Camera, (Pixel-Camera).normalize());
                pixelColor = pixelColor + scene.getColor(ray, max_depth);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., pow(pixelColor[0]/rays_per_pixel, 1./gamma)*255);
            image[(i * W + j) * 3 + 1] = std::min(255., pow(pixelColor[1]/rays_per_pixel, 1./gamma)*255);
            image[(i * W + j) * 3 + 2] = std::min(255., pow(pixelColor[2]/rays_per_pixel, 1./gamma)*255);
        }
    }
    // stop timer
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "Rendering took " << duration.count() << " milliseconds." << std::endl;

    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}