from django.db import models
from django.core.serializers.json import DjangoJSONEncoder
import numpy as np
import json

class NumpyEncoder(DjangoJSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)

class OperationStep(models.Model):
    session_id = models.PositiveIntegerField()
    step = models.CharField(max_length=100)
    operation_type = models.CharField(max_length=100)
    input_values = models.JSONField(encoder=NumpyEncoder)
    output_values = models.JSONField(null=True, blank=True, encoder=NumpyEncoder)

    def __str__(self):
        return f"Step {self.step}"