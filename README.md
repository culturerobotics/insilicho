![tests](https://github.com/culturerobotics/insilicho/actions/workflows/python_tests.yml/badge.svg)
![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)

# InSiliCHO
A model to capture CHO dynamics in-silico.  An accompanying streamlit app to try the model out is [available](https://culturebio-insilicho.streamlit.app/). 

Model is based on the following primary sources:
- Pörtner, Ralf, ed. Animal Cell Biotechnology: Methods and Protocols. Vol. 2095. Methods in Molecular Biology. New York, NY: Springer US, 2020. https://doi.org/10.1007/978-1-0716-0191-4. 
- Möller, Johannes, Tanja Hernández Rodríguez, Jan Müller, Lukas Arndt, Kim B. Kuchemüller, Björn Frahm, Regine Eibl, Dieter Eibl, and Ralf Pörtner. “Model Uncertainty-Based Evaluation of Process Strategies during Scale-up of Biopharmaceutical Processes.” Computers & Chemical Engineering 134 (March 2020): 106693. https://doi.org/10.1016/j.compchemeng.2019.106693.

Additional sources include:
- Parolini, Dott Nicola, and Susanna Carcano. “A model for cell growth in batch bioreactors,” 2009, Thesis.
- Frahm, Björn. “Seed Train Optimization for Cell Culture.” In Animal Cell Biotechnology, edited by Ralf Pörtner, 1104:355–67. Methods in Molecular Biology. Totowa, NJ: Humana Press, 2014. https://doi.org/10.1007/978-1-62703-733-4_22.

- Möller, Johannes, Kim B. Kuchemüller, Tobias Steinmetz, Kirsten S. Koopmann, and Ralf Pörtner. “Model-Assisted Design of Experiments as a Concept for Knowledge-Based Bioprocess Development.” Bioprocess and Biosystems Engineering 42, no. 5 (May 2019): 867–82. https://doi.org/10.1007/s00449-019-02089-7.




# Usage
This repo serves as a standalone package, available to install using (or added as a dependency) using:

`pip install insilicho`

# Example

```python

  from insilicho import run

  def T(time):
      """returns temperature in degC"""
      return 36.4

  def F(time):
      """returns flow rate in L/hr"""
      return 0.003

  model = run.GrowCHO(
      {"parameters": {"K_lys": "0.05 1/h"}},
      feed_fn=F,
      temp_fn=T,
  )

  model.execute(plot=True, initial_conditions={"V": "50 mL"})
  
  final_vol = model.full_result.state[-1, 8]
  print(final_V) # 0.914L

```
