defmodule AstroEquations.AstrophysicsAndAstronomy.InstrumentaionTest do
  use ExUnit.Case
  import AstroEquations.AstrophysicsAndAstronomy.Instrumentaion

  describe "lensmakers_equation/4-5" do
    test "calculates focal length with full equation" do
      assert_in_delta lensmakers_equation(1.5, 0.1, -0.1, 0.01), 0.10101, 1.0e-5
    end

    test "calculates focal length with thin lens approximation" do
      assert lensmakers_equation(1.5, 0.1, -0.1, 0.01, true) == 0.1
    end
  end

  describe "focal_ratio/2" do
    test "calculates f-number correctly" do
      assert focal_ratio(0.5, 0.1) == 5.0
    end
  end

  describe "field_of_view/3-4" do
    test "calculates FOV using f-number" do
      assert field_of_view(0.01, 0.2, 10) == 0.005
    end

    test "calculates FOV using system focal length" do
      assert field_of_view(0.01, 0.2, 10, 2.0) == 0.005
    end
  end

  describe "resolution limits" do
    @wavelength 500.0e-9
    @diameter 0.1
    @r0 0.2

    test "diffraction_limit/2" do
      assert diffraction_limit(@wavelength, @diameter) == 6.1e-6
    end

    test "seeing_limit/2" do
      assert seeing_limit(@wavelength, @r0) == 2.45e-6
    end

    test "total_resolution_limit/3" do
      result = total_resolution_limit(@wavelength, @diameter, @r0)
      assert_in_delta result, 6.58045e-6, 1.0e-8
    end
  end

  describe "nyquist_sampling/4" do
    test "calculates Nyquist sampling parameter" do
      assert nyquist_sampling(5.0e-6, 2.0, 500.0e-9, 0.2) == 10.0
    end
  end

  describe "plate_scale/1" do
    test "converts focal length to plate scale" do
      assert plate_scale(1.0) == {1.0, 206_265.0}
    end
  end

  describe "fitting_error/2" do
    test "calculates fitting error variance" do
      assert_in_delta fitting_error(0.1, 0.2), 0.102919, 1.0e-6
    end
  end

  describe "adaptive_optics_error/9" do
    test "calculates total adaptive optics error" do
      result = adaptive_optics_error(0.1, 0.2, 0.01, 0.02, 0.001, 0.002, 1.0, 500.0e-9, 10.0)
      assert_in_delta result, 0.11875, 1.0e-6
    end
  end

  describe "signal_to_noise/6" do
    test "calculates signal-to-noise ratio" do
      assert_in_delta signal_to_noise(100, 10, 5, 1000, 0.1, 2.0), 27.3861, 1.0e-4
    end
  end

  describe "atmospheric_extinction/3" do
    test "calculates apparent magnitude with extinction" do
      assert atmospheric_extinction(10.0, 0.2, :math.pi() / 3) == 10.4
    end
  end

  describe "tsiolkovsky_rocket_equation/3" do
    test "calculates delta-v using rocket equation" do
      assert_in_delta tsiolkovsky_rocket_equation(2500, 1000, 500), 1732.8679, 1.0e-4
    end
  end
end
