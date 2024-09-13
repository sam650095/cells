from django.db import models

class OperationStep(models.Model):
    session_id = models.PositiveIntegerField()
    step = models.CharField(max_length=100)
    operation_type = models.CharField(max_length=100)
    input_values = models.JSONField()
    output_values = models.JSONField(null=True, blank=True)

    def __str__(self):
        return f"Step {self.step}"