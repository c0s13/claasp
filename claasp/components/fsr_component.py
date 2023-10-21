
# ****************************************************************************
# Copyright 2023 Technology Innovation Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************


from sage.modules.free_module_element import vector

from claasp.input import Input
from claasp.component import Component, free_input


class FSR(Component):
    def __init__(self, current_round_number, current_round_number_of_components, input_id_links,
                 input_bit_positions, output_bit_size, description):
        component_id = f'fsr_{current_round_number}_{current_round_number_of_components}'
        component_type = 'fsr'
        input_len = 0
        for bits in input_bit_positions:
            input_len = input_len + len(bits)
        component_input = Input(input_len, input_id_links, input_bit_positions)
        super().__init__(component_id, component_type, component_input, output_bit_size, description)
        self.input_len = input_len

    def algebraic_polynomials(self, model):
        """
        Return a list of polynomials for the feedback shift registers.

        INPUT:

        - ``model`` -- **model object**; a model instance

        EXAMPLES::

            sage: from claasp.ciphers.stream_ciphers.a5_1_stream_cipher import A51StreamCipher
            sage: from claasp.cipher_modules.models.algebraic.algebraic_model import AlgebraicModel
            sage: a51 = A51StreamCipher()
            sage: fsr_component = a51.get_component_from_id("fsr_1_0")
            sage: algebraic = AlgebraicModel(a51)
            sage: L = fsr_component.algebraic_polynomials(algebraic)
            sage: L[0]
            linear_layer_0_6_y0 + linear_layer_0_6_x23 + linear_layer_0_6_x19 + linear_layer_0_6_x18 + linear_layer_0_6_x16 + linear_layer_0_6_x15 + linear_layer_0_6_x14 + linear_layer_0_6_x12 + linear_layer_0_6_x9 + linear_layer_0_6_x8 + linear_layer_0_6_x6 + linear_layer_0_6_x3
        """

        number_of_registers = self.description[0]
        output = BitArray(input)
        R = BooleanPolynomialRing(len(input), 'x')
        number_of_registers = len(registers_info)
        registers_polynomial = [0 for _ in range(number_of_registers)]
        registers_start = [0 for _ in range(number_of_registers)]
        registers_update_bit = [0 for _ in range(number_of_registers)]
        clock_polynomials = [None for _ in range(number_of_registers)]
        end = 0

        word_gf = GF(pow(2, bits_inside_word))
        word_array = bits_to_word(input, bits_inside_word, word_gf)
        R = PolynomialRing(word_gf, len(word_array), 'x')
        number_of_registers = len(registers_info)
        registers_polynomial = [0 for _ in range(number_of_registers)]
        registers_start = [0 for _ in range(number_of_registers)]
        registers_update_word = [0 for _ in range(number_of_registers)]
        clock_polynomials = [None for _ in range(number_of_registers)]
        end = 0

        noutputs = self.output_bit_size
        ninputs = self.input_bit_size
        ring_R = model.ring()
        x = list(ring_R, (map(ring_R, [self.id + "_" + model.input_postfix + str(i) for i in range(ninputs)])))
        y = vector(ring_R,
                   list(map(ring_R, [self.id + "_" + model.output_postfix + str(i) for i in range(noutputs)])))
        polynomial_index_list = self.description[0]
        loop = self.description[1]

        fsr_polynomial = 0
        for _ in polynomial_index_list:
            m = 1
            for __ in _ :
                m *= x[__]
            fsr_polynomial += m

        for i in range(loop):
            output_bit = fsr_polynomial(*x)
            x = x[:-1]
            x.append(output_bit)

        output_polynomials = y+vector(x)
        return output_polynomials


