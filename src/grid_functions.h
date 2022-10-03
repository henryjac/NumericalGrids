#include "domain.h"
#include "../../project2/src/matrix.h"

#include <memory>

typedef std::function<double(Point)> Fnc2D;

enum DiffCase {
    BOTTOM_LEFT, BOTTOM_RIGHT, TOP_LEFT, TOP_RIGHT,
    LEFT, RIGHT, BOTTOM, TOP,
    INSIDE
};

enum DiffDirection1D {
    XI, ETA
};

enum DiffDirection2D {
    XI_XI, ETA_ETA, XI_ETA
};

class GFkt : std::enable_shared_from_this<GFkt> {
    public:
        GFkt() : u(), grid(nullptr), h_eta(0), h_xi(0) {}
        GFkt(std::shared_ptr<Domain>);
        GFkt(std::shared_ptr<Domain>, Matrix);
        GFkt(const GFkt& U) : u(U.u), grid(U.grid), h_eta(U.h_eta), h_xi(U.h_xi) {}
        GFkt(GFkt&&);
        GFkt& operator=(const GFkt&);
        GFkt operator+(const GFkt&) const;
        GFkt operator-(const GFkt&) const;
        GFkt operator*(const GFkt&) const;
        GFkt operator*(double) const;
        GFkt operator/(double) const;

        void fill_matrix(Fnc2D);
        void printMatrix() {u.printMatrix();}
        void save(const char*, const char* gridfile=nullptr, const char* boundaryfile=nullptr);
        
        std::array<std::shared_ptr<GFkt>,2> pd();
        GFkt pdx();
        GFkt pdx2();
        GFkt pdy();
        GFkt pdy2();
        GFkt laplace();
        GFkt laplace2();
        
        double norm_1() {return u.norm_1();}
        double norm_inf() {return u.norm_inf();}
    private:
        Matrix u;
        std::shared_ptr<Domain> grid;
        double h_eta, h_xi;

        template <class intFnc>
        double pderiv(intFnc, DiffDirection1D, int, int);
        std::shared_ptr<Matrix> pd(int, int); 

        template <class intFnc>
        double pderiv2(intFnc, DiffDirection2D, int, int);
        std::shared_ptr<Matrix> pd2(int, int);

        DiffCase getDiffCase(int, int);
};
