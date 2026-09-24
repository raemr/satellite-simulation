#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

constexpr double pi = 3.14159265358979323846;
constexpr double earth_radius = 6371000.0;       // metres
constexpr double orbit_height = 550000.0;        // metres
constexpr double earth_mu = 3.986004418e14;      // m^3/s^2
constexpr double earth_spin = 7.2921159e-5;      // radians/s
constexpr double light_speed = 299792458.0;      // m/s
constexpr double carrier = 2.0e9;                // Hz
constexpr double beam_power = 10.0;               // W per beam
constexpr double peak_gain = 1000.0;              // linear gain
constexpr double beam_width = 2.0 * pi / 180.0;   // radians
constexpr double other_target_angle = 0.3 * pi / 180.0;
constexpr double bandwidth = 1.0e6;              // Hz
constexpr double noise_temperature = 290.0;      // K
constexpr double boltzmann = 1.380649e-23;       // J/K

struct Vec2 {
    double x;
    double y;
};

struct State {
    Vec2 position;
    Vec2 velocity;
};

Vec2 add(Vec2 a, Vec2 b) {
    return {a.x + b.x, a.y + b.y};
}

Vec2 subtract(Vec2 a, Vec2 b) {
    return {a.x - b.x, a.y - b.y};
}

Vec2 scale(Vec2 v, double factor) {
    return {factor * v.x, factor * v.y};
}

double dot(Vec2 a, Vec2 b) {
    return a.x * b.x + a.y * b.y;
}

double length(Vec2 v) {
    return std::sqrt(dot(v, v));
}

Vec2 unit(Vec2 v) {
    return scale(v, 1.0 / length(v));
}

State derivative(State state) {
    double radius = length(state.position);
    Vec2 acceleration = scale(state.position, -earth_mu / (radius * radius * radius));
    return {state.velocity, acceleration};
}

State add_scaled(State a, State b, double factor) {
    return {add(a.position, scale(b.position, factor)),
            add(a.velocity, scale(b.velocity, factor))};
}

State rk4_step(State state, double seconds) {
    State k1 = derivative(state);
    State k2 = derivative(add_scaled(state, k1, seconds / 2.0));
    State k3 = derivative(add_scaled(state, k2, seconds / 2.0));
    State k4 = derivative(add_scaled(state, k3, seconds));

    State next = add_scaled(state, k1, seconds / 6.0);
    next = add_scaled(next, k2, seconds / 3.0);
    next = add_scaled(next, k3, seconds / 3.0);
    return add_scaled(next, k4, seconds / 6.0);
}

Vec2 target_position(double time, double starting_angle) {
    double angle = starting_angle + earth_spin * time;
    return {earth_radius * std::cos(angle), earth_radius * std::sin(angle)};
}

Vec2 target_velocity(double time, double starting_angle) {
    Vec2 position = target_position(time, starting_angle);
    return {-earth_spin * position.y, earth_spin * position.x};
}

double travel_time(State satellite, double send_time, double starting_angle) {
    double delay = length(subtract(target_position(send_time, starting_angle), satellite.position))
                   / light_speed;
    for (int i = 0; i < 4; ++i) {
        delay = length(subtract(target_position(send_time + delay, starting_angle),
                                satellite.position)) / light_speed;
    }
    return delay;
}

void show_link(State satellite, double send_time) {
    double delay = travel_time(satellite, send_time, 0.0);
    Vec2 target = target_position(send_time + delay, 0.0);
    Vec2 old_target = target_position(send_time, 0.0);
    double no_lead_delay = length(subtract(old_target, satellite.position)) /
                           light_speed;
    double no_lead_miss = length(subtract(target_position(send_time + no_lead_delay, 0.0),
                                            old_target));
    Vec2 toward_target = unit(subtract(target, satellite.position));
    Vec2 toward_satellite = scale(toward_target, -1.0);
    double distance = length(subtract(target, satellite.position));

    double sine_elevation = dot(unit(target), toward_satellite);
    double elevation = std::asin(std::clamp(sine_elevation, -1.0, 1.0)) * 180.0 / pi;

    double ground_speed_along_ray = dot(toward_target,
                                       target_velocity(send_time + delay, 0.0));
    double satellite_speed_along_ray = dot(toward_target, satellite.velocity);
    double frequency_ratio = (light_speed - ground_speed_along_ray) /
                             (light_speed - satellite_speed_along_ray);
    double received_shift = carrier * (frequency_ratio - 1.0);
    double transmit_shift = carrier * (1.0 / frequency_ratio - 1.0);

    double other_delay = travel_time(satellite, send_time, other_target_angle);
    Vec2 other_target = target_position(send_time + other_delay, other_target_angle);
    Vec2 other_beam = unit(subtract(other_target, satellite.position));
    double angle_between_beams = std::acos(std::clamp(dot(toward_target, other_beam),
                                                     -1.0, 1.0));
    double other_gain = peak_gain * std::exp(-0.5 * std::pow(angle_between_beams /
                                                              beam_width, 2.0));
    double wavelength = light_speed / carrier;
    double path_factor = std::pow(wavelength / (4.0 * pi * distance), 2.0);
    double desired_power = beam_power * peak_gain * path_factor;
    double interfering_power = beam_power * other_gain * path_factor;
    double noise_power = boltzmann * noise_temperature * bandwidth;
    double sinr_db = 10.0 * std::log10(desired_power /
                                      (interfering_power + noise_power));

    bool visible = elevation >= 10.0;
    std::cout << send_time << "  " << distance / 1000.0 << "  "
              << elevation << "  " << delay * 1000.0 << "  " << no_lead_miss << "  "
              << received_shift << "  " << transmit_shift << "  ";
    if (visible) {
        std::cout << sinr_db;
    } else {
        std::cout << "--";
    }
    std::cout << "  " << (visible ? "yes" : "no") << '\n';
}

int main() {
    double orbit_radius = earth_radius + orbit_height;
    double orbital_speed = std::sqrt(earth_mu / orbit_radius);
    double angular_speed = orbital_speed / orbit_radius;
    double start_time = -360.0;
    double start_angle = angular_speed * start_time;

    State satellite = {
        {orbit_radius * std::cos(start_angle), orbit_radius * std::sin(start_angle)},
        {-orbital_speed * std::sin(start_angle), orbital_speed * std::cos(start_angle)}
    };

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "time_s  range_km  elev_deg  travel_ms  no_lead_m  doppler_Hz  tx_fix_Hz  sinr_dB  visible\n";
    for (int step = 0; step <= 720; ++step) {
        double time = start_time + step;
        if (step % 120 == 0) {
            show_link(satellite, time);
        }
        if (step < 720) {
            satellite = rk4_step(satellite, 1.0);
        }
    }

    Vec2 expected = {orbit_radius * std::cos(angular_speed * 360.0),
                     orbit_radius * std::sin(angular_speed * 360.0)};
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "orbit check: final position differs from a circle by "
              << length(subtract(satellite.position, expected)) << " m\n";
}
