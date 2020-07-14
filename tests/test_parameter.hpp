using namespace cvx;

TEST_CASE("Parameter")
{
    double zero = 0.;
    double one = 1.;
    double two = 2.;
    double three = 3.;

    { // Constants
        Parameter p0(zero);
        Parameter p1(one);
        Parameter p2(two);
        Parameter p3(three);

        REQUIRE(p0 == p0);
        REQUIRE(p1 == p1);
        REQUIRE(p2 == p2);
        REQUIRE_FALSE(p0 == p1);
        REQUIRE_FALSE(p2 == p0);
        REQUIRE_FALSE(p2 == p1);

        Parameter p0p1 = p0 + p1;
        Parameter p1p2 = p1 + p2;
        Parameter p2p3 = p2 + p3;

        Parameter p0t1 = p0 * p1;
        Parameter p1t2 = p1 * p2;
        Parameter p2t3 = p2 * p3;

        // Addition
        REQUIRE((p0 + p1).getValue() == 1.);
        REQUIRE((p0 + p2).getValue() == 2.);
        REQUIRE((p1 + p2).getValue() == 3.);

        // Multiplication
        REQUIRE((p0 * p1).getValue() == 0.);
        REQUIRE((p1 * p0).getValue() == 0.);
        REQUIRE((p0 * p2).getValue() == 0.);
        REQUIRE((p2 * p0).getValue() == 0.);
        REQUIRE((p3 * p1).getValue() == 3.);
        REQUIRE((p1 * p3).getValue() == 3.);
        REQUIRE((p2 * p3).getValue() == 6.);
        REQUIRE((p1 * p3).getValue() == 3.);

        // Square root
        REQUIRE((sqrt(p0).getValue()) == 0.);
        REQUIRE((sqrt(p1).getValue()) == 1.);
        REQUIRE((sqrt(p2).getValue()) == std::sqrt(2.));
        REQUIRE((sqrt(p2p3).getValue()) == std::sqrt(5.));
        REQUIRE((sqrt(p2t3).getValue()) == std::sqrt(6.));
    }

    { // Pointers
        Parameter p0(&zero);
        Parameter p1(&one);
        Parameter p2(&two);
        Parameter p3(&three);

        REQUIRE(p0 == p0);
        REQUIRE(p1 == p1);
        REQUIRE(p2 == p2);
        REQUIRE_FALSE(p0 == p1);
        REQUIRE_FALSE(p2 == p0);
        REQUIRE_FALSE(p2 == p1);

        Parameter p0p1 = p0 + p1;
        Parameter p1p2 = p1 + p2;
        Parameter p2p3 = p2 + p3;

        Parameter p0t1 = p0 * p1;
        Parameter p1t2 = p1 * p2;
        Parameter p2t3 = p2 * p3;

        // Addition
        REQUIRE((p0 + p1).getValue() == 1.);
        REQUIRE((p0 + p2).getValue() == 2.);
        REQUIRE((p1 + p2).getValue() == 3.);

        // Multiplication
        REQUIRE((p0 * p1).getValue() == 0.);
        REQUIRE((p1 * p0).getValue() == 0.);
        REQUIRE((p0 * p2).getValue() == 0.);
        REQUIRE((p2 * p0).getValue() == 0.);
        REQUIRE((p3 * p1).getValue() == 3.);
        REQUIRE((p1 * p3).getValue() == 3.);
        REQUIRE((p2 * p3).getValue() == 6.);
        REQUIRE((p1 * p3).getValue() == 3.);

        // Square root
        REQUIRE((sqrt(p0).getValue()) == 0.);
        REQUIRE((sqrt(p1).getValue()) == 1.);
        REQUIRE((sqrt(p2).getValue()) == std::sqrt(2.));
        REQUIRE((sqrt(p2p3).getValue()) == std::sqrt(5.));
        REQUIRE((sqrt(p2t3).getValue()) == std::sqrt(6.));

        REQUIRE(p0 * p0 == p0 * p0);
        REQUIRE(p1 + p0 == p0 + p1);
        REQUIRE(p2 * p2 == p2 * p2);
        REQUIRE_FALSE(p0 * p1 == p1 * p2);
        REQUIRE_FALSE(p2 + p1 == p2 * p1);
        REQUIRE_FALSE(p2 == p1 * p2);
    }
}
